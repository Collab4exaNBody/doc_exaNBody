ExaNBody: Framework for N-Body Simulations on HPC Platforms
===========================================================

This framework, developed at the French Atomic Agency (CEA), is tailored for N-Body simulations on High-Performance Computing (HPC) platforms. Originally designed for the ExaSTAMP Molecular Dynamics code, it has been extended to cater to various N-Body problems.

Key Characteristics:

* **Language:** Implemented in C++17, leveraging modern language features for efficiency and versatility.

* **Parallelization:**: Hybrid approach integrating:
	* Vectorization for CPU optimization.
	* Thread-parallelization using OpenMP for multi-core architectures.
	* GPU-parallelization via CUDA to harness GPU computational power.
	* MPI-parallelization for distributed memory systems.
  
* **Spatial Domain Decomposition:** Utilizes spatial domain decomposition techniques for efficient workload distribution among processors.

* **Load Balancing (RCB):** Implements Load Balancing using Recursive Coordinate Bisection (RCB) for optimal task distribution among processing units.

* **Parallel IO:** Enables efficient handling of parallel Input/Output operations. Supports checkpoint files, parallel Paraview files, and diagnostics in a parallelized manner.

* **In-situ Analysis:** Provides real-time data analysis capabilities during simulation execution, minimizing data movement and storage overhead.


ExaDEM Variant
--------------

ExaSPH Variant
--------------

ExaSTAMP Variant
----------------
