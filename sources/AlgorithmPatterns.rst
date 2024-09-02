
Algorithm Patterns
==================

This is Developer Documentation.


The aim of this section is to explain to developers the tools used in exaNBody to propose default parallel algorithms. For this purpose, we will describe the 4 main algorithmic patterns:

- `block_parallel_for`: exposes a simple parallelization scheme (like a parallel for in OpenMP), but with a notion of "thread blocks" and the possibility of asynchronous execution. It is defined directly in the Onika layer, as it is not specific to particle systems.
- `compute_cell_particles`: performs a parallel processing on all particles contained within a grid (Grid type pattern).
- `reduce_cell_particles`: performs an operation on each particle returning a value, the values obtained being reduced to a single one (e.g., a sum).
- `compute_cell_particle_pairs`: the most important and most used pattern, allows calculations on particles and their neighbors.

Guidelines:

- Support hybrid execution, automatically directing calculations either to the `CPU` cores or to the `GPU`.
- Avoid code duplication between components.
- Avoid the complexities of setting up `OpenMP`/`CUDA` parallelism.
- Allow asynchronous parallel operations (with results retrieved later in the code).
- Automatically adapt to the different `OpenMP` parallelism modes supported by the component system (simple, nested, and task-based).
- Facilitate the implementation of a single computation function for both `CPU` and `GPU`.
- Manage dynamic load balancing among the `SMs` of the `GPUs` themselves.

Block Parallel For
------------------

This construction is the most rudimentary; it is equivalent to an `OpenMP` parallel for directive but with many additional possibilities. We will look at a step-by-step example where a fixed value is added to all elements of an array of arrays. The example on a two-dimensional array demonstrates the use of two-level multi-threaded parallelism we have discussed, which resembles a bit the array of cells, each containing a set of particles typically processed in exaNBody.


Basic usage
^^^^^^^^^^^

Before calling `block_parallel_for`, it is necessary to define the functor (the calculation function) that will be applied in parallel to our data.

.. code-block:: cpp

  namespace exaStamp
  {
    using namespace exanb;
    // a functor = function applied in parallel
    // it is defined by a class (or struct) with the call operator ()
    struct ValueAddFunctor
    {
      double **array;      // our 2D array
      size_t columns;      // number of columns in our array
      double m_value_to_add; // value to add
     
      ONIKA_HOST_DEVICE_FUNC       // works on CPU and GPU
      void operator () (size_t i)  // call operator with i in [0;n[
      {                            // a whole block (all its threads) execute iteration i
        ONIKA_CU_BLOCK_SIMD_FOR(size_t, j, 0, columns)  // parallelization among the threads of the current block
        {                                               // for iterations on j in [0;columns[
          array[i][j] += m_value_to_add;               // each thread executes 0, 1, or multiple iterations of j
        }
      }
    };
  }


- `ONIKA_CU_BLOCK_SIMD_FOR`: Defines a macro that facilitates the use of the second level of parallelism within threads of the same `CUDA` block.
- Usage: 

  - GPU :Ensures that the functor's code is executed across all threads within a `CUDA` block with consistent iteration index `i`. When compiled for the `GPU`, it distributes iterations based on `ONIKA_CU_BLOCK_SIZE` and `ONIKA_CU_THREAD_IDX`. 
  - CPU: On the host, this construct transforms into a straightforward for loop. However, on the host, the loop includes the `OpenMP` directive #pragma omp simd for, ensuring loop vectorization, even with only one  `CPU` thread per block, provided `CPU` and compiler support is available.

Next, we need to inform the compiler whether our functor is compatible with `CUDA`. In some situations, our code may have a strong incompatibility (such as calling a non-compatible library) that would prevent it from compiling correctly with the nvcc layer. This needs to be known at compile time (not just at runtime) because otherwise, compilation errors would occur.

.. code-block:: cpp

  namespace onika
  {
    namespace parallel
    {
      template<>
      struct BlockParallelForFunctorTraits<ValueAddFunctor>
      {
        static constexpr bool CudaCompatible = true; // or false to prevent the code from being compiled with CUDA
      };
    }
  }

The final step is to launch the parallel operation in the code of our component:

.. code-block:: cpp

  namespace exaStamp 
  {
    class MyValueAdd : public OperatorNode
    {
      ADD_SLOT(Array2D, my_array, INPUT, REQUIRED);
      ADD_SLOT(double, my_value, INPUT, 1.0);
      public:
      inline void execute() override final {
        ValueAddFunctor func = { my_array->data()     // data of our array (double**), initializes the 'array' member of our functor
                               , my_array->columns()  // number of columns in our 2D array
                               , *my_value };         // value to add to the elements of the array
        // Launching the parallel operation, which can execute on GPU if the execution context allows
        onika::parallel::block_parallel_for(my_array->rows()             // number of iterations, parallelize at the first level over rows
                                           , func                         // the function to call in parallel
                                           , parallel_execution_context() // returns the parallel execution context associated with this component
                                           );
      }
    };
  }



Asynchrone usage
^^^^^^^^^^^^^^^^

In the following example, we are still using the simplest version of `block_parallel_for`. Optional parameters can be added to control asynchronicity or even disable execution on the `GPU`.

.. code-block:: cpp

  namespace exaStamp 
  {
    class MyValueAdd : public OperatorNode
    {
      ADD_SLOT(Array2D, my_array, INPUT, REQUIRED);
      ADD_SLOT(double, my_value, INPUT, 1.0);
  public:
      inline void execute() override final {
        size_t cols = my_array->columns();
        size_t rows = my_array->rows();
        bool enable_gpu = (cols * rows) > 1000;     // enable GPU execution only if the matrix has more than 1000 values
        ValueAddFunctor func = { my_array->data(), cols, *my_value };
        auto control = onika::parallel::block_parallel_for(
                         rows
                       , func
                       , parallel_execution_context()
                       , enable_gpu                    // enable or disable GPU execution
                       , true                          // request asynchronous execution
                       );
        std::cout << "Meanwhile, the parallel operation is executing..." << std::endl;
        control->wait();                               // wait for the operation to complete and results to be ready to read
        std::cout << "Parallel operation completed!" << std::endl;
      }
    };
  }


If the execution context allows it, the parallel operation will proceed in the background, occupying either the free threads (other than those executing this code) or the `GPU`. This can be very useful, especially for overlapping computations and `MPI` message sends. 

.. warning::

  The operation is asynchronous only if the execution context permits it. Otherwise, it will proceed synchronously and complete before the `block_parallel_for` function returns. In such cases, calling `control->wait()` will simply have no effect.


Compute Cell Particles
----------------------

This is the first specialized algorithmic pattern for particle systems. Therefore, it can only be applied to an instantiation of the `Grid` type from `exaNBody`. This pattern is the simplest; it will perform an independent operation in parallel on all particles in the system. To understand how this works in the implementation of a component, the following code presents a fully commented example aimed at increasing the velocity of all particles by a constant vector.

.. code-block:: cpp

  #include <exanb/core/operator.h>
  #include <exanb/core/operator_slot.h>
  #include <exanb/core/operator_factory.h>
  #include <exanb/core/make_grid_variant_operator.h> // because we use make_grid_variant_operator for component registration
  #include <exanb/compute/compute_cell_particles.h>  // provides the compute_cell_particles function

  namespace exaStamp {
    struct AddVec3Functor        // our functor, adds a vector to the 3 components of particle velocity
    {
      const Vec3d vec_to_add;    // functor parameter: the velocity to add

      ONIKA_HOST_DEVICE_FUNC
      inline void operator ()(double& vx, double& vy, double& vz) const // parameters defined by compute_fields_t declared below
      {
        vx += vec_to_add.x;  // perform the desired operation
        vy += vec_to_add.y;  // particle attributes are passed by reference (&), so they can be modified
        vz += vec_to_add.z;  // here, we only need to implement the local operation on the particle, without worrying about thread block levels
      }
    };
  }

  namespace exanb {
    template<> struct ComputeCellParticlesTraits<exaStamp::AddVec3Functor>
    {
      static inline constexpr bool RequiresBlockSynchronousCall = false; // no collaboration between threads of the same block
      static inline constexpr bool CudaCompatible = true;                // compatible with Cuda (thanks to ONIKA_HOST_DEVICE_FUNC usage)
    };
  }

  namespace exaStamp 
  {
    template< class GridT                                                        // our component adapts to any grid type
            , class = AssertGridHasFields<GridT,field::_vx,field::_vy,field::_vz> > // as long as particles have the vx, vy, and vz attributes
    class ParticleVelocityAdd : public OperatorNode
    {
        ADD_SLOT(Vec3d, vec, INPUT, REQUIRED);    // vector to add to particle velocities
        ADD_SLOT(GridT, grid, INPUT_OUTPUT);      // grid containing cells and all particles within the sub-domain
  public:
      inline void execute() override final {
        using compute_fields_t = FieldSet<field::vx, field::vy, field::vz>; // fields on which our functor operates
        AddVec3Functor func = { *vec };                                     // instantiate the functor with the velocity provided as input to the component
        compute_cell_particles(*grid                    // grid containing particles
                              , false                    // do not apply our function in ghost zones
                              , func                     // the functor to apply
                              , compute_fields_t{}       // attributes on which to compute ⇒ defines the functor's call parameters
                              , parallel_execution_context() // component's parallel execution context
                              );
      }
    };
    // component registration
    template<class GridT> using ParticleVelocityAddTmpl = ParticleVelocityAdd<GridT>;
    CONSTRUCTOR_FUNCTION {
      OperatorNodeFactory::instance()->register_factory("add_velocity", make_grid_variant_operator<ParticleVelocityAddTmpl>);
    }
  }

.. note::

	To improve the performance of Compute Cell Particles, you can choose to run it only on non-empty cells, meaning cells that contain at least one particle. This feature is particularly useful in DEM (Discrete Element Method), where it is common 		to encounter a significant number of empty cells. To use this feature, you can rely on the default parameters `filled_cells`, which is a list containing the indexes of the non-empty cells, and `number_filled_cells`, which indicates the number of non-empty cells.


Reduce Cell Particles
---------------------

This algorithm pattern is very similar to the previous one, but it allows for computing the reduction of a calculation result across all particles. As we will see in the commented example below, the usage of this construction will be very similar to `compute_cell_particles` with one important exception: the functor can be called in multiple ways (multiple implementations of the operator () are possible). This reflects the need to perform reduction in three stages for performance reasons: local reduction within each thread, reduction of partial sums computed by threads within the same thread block, and finally the overall reduction that sums the partial contributions from different thread blocks.

.. code-block:: cpp

  #include <exanb/core/operator.h>
  #include <exanb/core/operator_slot.h>
  #include <exanb/core/operator_factory.h>
  #include <exanb/core/make_grid_variant_operator.h> // because we use make_grid_variant_operator for component registration
  #include <exanb/compute/reduce_cell_particles.h>   // provides the reduce_cell_particles function

  namespace exaStamp {
    struct ReduceVec3NormFunctor // our functor calculates the sum of norms of forces
    {
      ONIKA_HOST_DEVICE_FUNC                                    // operator for local reduction within a thread
      inline void operator ()( double & sum_norm                // reference to accumulate contributions
                             , double fx, double fy, double fz  // particle force, parameters determined by reduce_fields_t declared below
                             , reduce_thread_local_t            // phantom parameter to differentiate call forms (here local thread reduction)
                             ) const {
        sum_norm += sqrt( fx*fx + fy*fy + fz*fz );              // compute norm and add contribution
      }

      ONIKA_HOST_DEVICE_FUNC                          // operator for internal reduction within a thread block to a single block value
      inline void operator ()( double& sum_norm       // reference to accumulate contributions from threads within the block
                             , double other_sum       // one of the partial sums to accumulate
                             , reduce_thread_block_t  // indicates block-level reduction
                             ) const {
        ONIKA_CU_ATOMIC_ADD( sum_norm , other_sum );  // atomic addition function (thread-safe), works in both CUDA and CPU context
      }

      ONIKA_HOST_DEVICE_FUNC                       // operator for reduction across all thread blocks
      inline void operator ()( double& sum_norm    // reference to global result
                             , double other_sum    // contribution from one block to add
                             , reduce_global_t     // indicates final reduction
                             ) const {
        ONIKA_CU_ATOMIC_ADD( sum_norm , other_sum );
      }
    };
  }

  namespace exanb {
    template<> struct ReduceCellParticlesTraits<exaStamp::ReduceVec3NormFunctor>
    {
      static inline constexpr bool RequiresBlockSynchronousCall = false; // does not use intra-block thread collaboration
      static inline constexpr bool RequiresCellParticleIndex = false;    // no additional parameters (cell/particle indices)
      static inline constexpr bool CudaCompatible = true;                // CUDA compatible
    };
  }

  namespace exaStamp {
  template< class GridT                                                      // our component adapts to any grid type
          , class = AssertGridHasFields<GridT,field::_fx,field::_fy,field::_fz> > // as long as particles have fx, fy, and fz attributes
  class SumForceNorm : public OperatorNode
  {
      ADD_SLOT(GridT, grid, INPUT, REQUIRED);    // grid containing cells and all particles within the sub-domain
      ADD_SLOT(double, sum_norm, OUTPUT);       // output value = sum of force norms
  public:
    inline void execute() override final {
      using reduce_fields_t = FieldSet<field::fx, field::fy, field::fz>; // fields on which our functor operates
      ReduceVec3NormFunctor func = {};                                    // instantiate the functor for summing force norms
      *sum_norm = 0.0;
      reduce_cell_particles(*grid                         // grid containing particles of the system
                           , false                         // do not compute in ghost zones
                           , func                          // functor to use
                           , *sum_norm                     // initial value for reduction input, final result output
                           , reduce_fields_t{}             // fields used for reduction, defines call parameters
                           , parallel_execution_context()  // current component's parallel execution context
                           );
    }
  };
  
  // component registration similar to the one in the example for compute_cell_particles
  }


This code defines a component (``SumForceNorm``) that computes the sum of norms of forces acting on particles within a grid (`grid`). Key points to note:

- It uses a functor (ReduceVec3NormFunctor) with multiple operator overloads to perform reduction across particles in three stages: local to each thread, within each thread block, and globally across all thread blocks.
- Traits (``ReduceCellParticlesTraits``) are specialized to specify that the functor supports `CUDA` and does not require intra-block thread synchronization.
- The component is registered using ``OperatorNodeFactory``, making it accessible under a specific name for instantiation.

This example demonstrates advanced parallel computation techniques within a particle system framework (`exaNBody`), focusing on efficient reduction operations across large datasets.

.. note::

	To improve the performance of Compute Cell Particles, you can choose to run it only on non-empty cells, meaning cells that contain at least one particle. This feature is particularly useful in DEM (Discrete Element Method), where it is common 		to encounter a significant number of empty cells. To use this feature, you can rely on the default parameters `cells_idx`, which is a list containing the indexes of the non-empty cells, and `n_cells`, which indicates the number of non-empty cells.


Compute Pair Interaction
------------------------

Pattern Overview:

- Purpose: Compute interactions (potentials or forces) between pairs of particles based on their proximity (rij < rcut).
  - Components: Utilizes both particle grids and neighbor lists (see neighbor lists section).
  - Two Invocation Modes:

    - With Buffer: Collects attributes of neighboring particles into a buffer before invoking the user-defined functor once with this buffer.
    - Without Buffer: Computes interactions "on-the-fly" without pre-collecting neighbors' attributes.

- Usage Examples:
  - Scenario 1: Interaction potentials require knowledge of all neighboring particles of at least one of the particles involved. Here, `compute_cell_particle_pairs` first accumulates neighbors' attributes into a buffer and then invokes the user functor.
  - Scenario 2: Interaction potentials only require attributes of the two particles involved, making it efficient to compute interactions without using an intermediate buffer.

- Flexibility:
  - Functors can be implemented to support both invocation protocols (with and without buffer), allowing the subsystem to choose the optimal method based on execution environment and computational requirements.
  - Functors can choose to implement only one of these protocols if the scenario allows for a straightforward implementation.

- Implementation:
  - The pattern involves defining a functor (``ComputePairInteractionFunctor``) that encapsulates the logic for computing interactions between pairs of particles.
  - Traits (``ComputeCellParticlePairsTraits``) are used to specify compatibility with `CUDA` and whether synchronous block calls are required.


Header files:

.. code-block:: cpp

  #pragma xstamp_enable_cuda                           // Enable compilation with nvcc, allowing GPU code generation
  #include <exanb/core/operator.h>                     // Base to create a component
  #include <exanb/core/operator_factory.h>             // For registering the component
  #include <exanb/core/operator_slot.h>                // For declaring component slots
  #include <exanb/core/make_grid_variant_operator.h>   // For registering a template component

  #include <exanb/core/grid.h>                         // Defines the Grid type containing particles
  #include <exanb/core/domain.h>                       // Domain type, represents the simulation domain 
  #include <exanb/particle_neighbors/chunk_neighbors.h>// GridChunkNeighbors type ⇒ lists of neighbors
  #include <exanb/compute/compute_cell_particle_pairs.h>// For compute_cell_particle_pairs function


Functor Creation:

For this example (lennard jones potential), we implement both possible forms of invocation so that our functor can be used with or without an intermediate buffer. This allows `compute_cell_particle_pairs` to choose the most efficient method based on the context. In practice, you can implement only one of these forms (if your potential formulation allows) or both-either way, it will work on both `CPU` and `GPU`.

.. code-block:: cpp

  namespace microStamp {
    using namespace exanb;
    
    struct LennardJonesForceFunctor {
      const LennardJonesParms m_params;  // parameters of our interaction potential
      // Operator for using an intermediate buffer containing all neighbors
      template<class ComputePairBufferT, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC
      inline void operator () (
          size_t n,                                // number of neighboring particles
          const ComputePairBufferT& buffer,         // intermediate buffer containing neighbor information
          double& e,                                // reference to energy, where we accumulate energy contribution
          double& fx, double& fy, double& fz,       // references to force components where we add interaction contribution
          CellParticlesT* cells                    // array of all cells, in case additional information is needed
      ) const {
        double _e = 0.0;                           // local contributions, initialized to 0
        double _fx = 0.0;
        double _fy = 0.0;
        double _fz = 0.0;
      
        for (size_t i = 0; i < n; ++i) {           // loop over neighbors in the buffer
          const double r = std::sqrt(buffer.d2[i]);
          double pair_e = 0.0, pair_de = 0.0;
        
          lj_compute_energy(m_params, r, pair_e, pair_de);  // calculate energy and its derivative
        
          const auto interaction_weight = buffer.nbh_data.get(i);
          pair_e *= interaction_weight;                    // weight the interaction and normalize over distance rij
          pair_de *= interaction_weight / r;              
        
          _fx += pair_de * buffer.drx[i];                  // add contributions from the i-th neighbor
          _fy += pair_de * buffer.dry[i];
          _fz += pair_de * buffer.drz[i];
          _e += 0.5 * pair_e;
        }
      
        e += _e;   // add local contributions to energy and force fields of the central particle
        fx += _fx;
        fy += _fy;
        fz += _fz;
      }

      // Operator for without buffer: one call per neighboring particle
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC
      inline void operator () (
          Vec3d dr,                           // relative position of the neighboring particle
          double d2,                          // square of the distance to the neighboring particle
          double& e, double& fx, double& fy, double& fz,  // references to variables to update with energy and force contributions of the pair
          CellParticlesT* cells,              // array of all cells, in case additional information is needed
          size_t neighbor_cell,               // index of the cell where the neighboring particle resides
          size_t neighbor_particle,           // index of the neighboring particle within its cell
          double interaction_weight           // weighting to apply on the interaction
      ) const {
        const double r = sqrt(d2);
        double pair_e = 0.0, pair_de = 0.0;
      
        lj_compute_energy(m_params, r, pair_e, pair_de);  // calculate energy and its derivative
        pair_e *= interaction_weight;                    // weight and normalize by distance
        pair_de *= interaction_weight / r;              
        
        fx += pair_de * dr.x;                            // add contributions of the pair to energy and forces of the central particle
        fy += pair_de * dr.y;
        fz += pair_de * dr.z;
        e += 0.5 * pair_e;
      }
    };
  }

Define the compile-time characteristics of the functor:

.. code-block:: cpp

  namespace exanb
  {
    template<> struct ComputePairTraits< microStamp::LennardJonesForceFunctor > // specialization for our functor
    {
      static inline constexpr bool RequiresBlockSynchronousCall = false; // no collaboration between threads within a block
      static inline constexpr bool ComputeBufferCompatible = true;  // allows invocation with a compute buffer
      static inline constexpr bool BufferLessCompatible = true;     // allows invocation without a buffer, for each neighbor
      static inline constexpr bool CudaCompatible = true;           // compatible with Cuda
    };
  }


Performing the calculation in parallel within the execute method of our component:

.. code-block:: cpp

  template<class GridT,                                                            // our component adapts to any type of grid
    class=AssertGridHasFields<GridT,field::_ep,field::_fx,field::_fy,field::_fz> > // which particles have attributes ep, fx, fy, and fz
  class LennardJonesForce : public OperatorNode {
    ADD_SLOT( LennardJonesParms        , config         , INPUT, REQUIRED);  // parameters of our potential
    ADD_SLOT( double                   , rcut           , INPUT, 0.0);       // cutoff radius
    ADD_SLOT( exanb::GridChunkNeighbors, chunk_neighbors, INPUT, exanb::GridChunkNeighbors{}); // neighbor lists
    ADD_SLOT( bool                     , ghost  , INPUT , false);            // indicates if we calculate in ghost zones
    ADD_SLOT( Domain                   , domain , INPUT , REQUIRED);         // simulation domain
    ADD_SLOT( GridT                    , grid   , INPUT_OUTPUT );            // grid containing particles in the local subdomain
    using ComputeBuffer = ComputePairBuffer2<>;                              // shortcut for the type "buffer containing neighbors"
    using CellParticles = typename GridT::CellParticles;                     // shortcut for the type "grid cell"
    using NeighborIterator = exanb::GridChunkNeighborsLightWeightIt<false>;  // shortcut for the type "neighbor iterator"
    using ComputeFields = FieldSet<field::_ep    // type defining the list of fields on which our functor operates
                                  ,field::_fx    // this list determines part of the functor's call parameters
                                  ,field::_fy    // we declare ep, fx, fy, and fz, which are of type double, so we will find
                                  ,field::_fz>;  // references to these attributes as parameters of our functor.
  public:
    inline void execute () override final
    {
      auto optional = make_compute_pair_optional_args(        // builds the list of particle traversal options
                        NeighborIterator{ *chunk_neighbors }  // iterator over neighbor lists
                      , ComputePairNullWeightIterator{}       // iterator over weighting values (no weighting here)
                      , LinearXForm{ domain->xform() }        // transformation to apply on particle positions (domain deformation here)
                      , ComputePairOptionalLocks<false>{});   // accessor to particle access locks (no locking here)
      compute_cell_particle_pairs( *grid                                     // grid containing particles
                                 , *rcut                                     // cutoff radius
                                 , *ghost                                    // whether to calculate in ghost zones or not
                                 , optional                                  // particle traversal options
                                 , make_compute_pair_buffer<ComputeBuffer>() // creates neighbor buffers
                                 , LennardJonesForceFunctor{ *config }       // instantiation of our functor with user parameters
                                 , ComputeFields{}                           // attributes on which our functor operates
                                 , DefaultPositionFields{}                   // uses default attributes (rx, ry, and rz) for particle positions
                                 , parallel_execution_context() );           // parallel execution context of the current component
    }
  };
