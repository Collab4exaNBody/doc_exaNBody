ExaNBody Command Lines
======================

ExaNBody offers a range of command-line options to facilitate debugging, profiling, and configuring OpenMP settings, providing users with enhanced control and flexibility over the behavior and performance of the code. These options enable developers to diagnose issues, analyze performance characteristics, and optimize parallel execution on multicore architectures. From enabling debug mode to specifying thread affinity and scheduling policies, users can tailor ExaNBody's behavior to suit their specific requirements and achieve optimal performance. This section outlines the various command-line options available, empowering users to leverage ExaNBody's capabilities effectively.

Command line and input file interaction
---------------------------------------

exaNBody based apps treat command line the same way as input files, the command line just being YAML elements expressed with another syntax.
YAML docuement, as processed by an exaNBody app is formed by a set of included files (either implicitly or explcitly through the 'includes' list) and the user input file passed
as the first argument of command line, as in the following exemple :
.. code-block::bash
  ./exaStamp myinput.msp

The YAML docuement built from user input and its included files have 3 reserved dictionary entries, namely 'configuration', reserved for configuring execution sub-system, 'includes', reserved to include other YAML files
and 'simulation' wich is interpreted as the the root batch operator representing the whole simulation to run.

When it comes to interpreting command line arguments, the exaNBody based application processes it abuild a YAML document out of it, merging it with the user defined YAML docuement as if it would have been included after all other included files.
Command line arguments are converted to a YAML document using the following conventions :

* If argument starts with '--set-', it is understood as a generic YAML dictionary entry, and each '-' is interpreted as a marker for an inner dictionary entry. The following value is understood as beeing the value associated with the given dictionary key. If no value found after a '--set-xxx' style argument, the value 'true' is implicitly used. For instance, --set-global-rcut_inc '1.2 m' is equivalent to including a YAML file wich contains
.. code-block::yaml

  global:
    rcut_inc: 1.2 m

* If a command line argument starts with --xxx with xxx being anything but 'set' , then similar rules apply as for the --set-xxx args, but the YAML block is understood as being a dictionary entry inside the 'configuration' block. For instance, --logging-debug true is equivalent to inclusion of a YAML file containing :
.. code-block::yaml

  configuration:
  logging:
    debug: true


Tune your run with OpenMP
-------------------------

Harnessing the power of OpenMP parallelization, ExaNBody provides users with the ability to fine-tune their execution environment for optimal performance on multicore architectures. Through a selection of command-line options, users can customize thread management, affinity settings, and loop scheduling to maximize parallel efficiency. This subsection introduces the command-line options available for configuring OpenMP behavior within ExaNBody.


.. list-table:: ExaNBody OpenMP Command Lines 
  :widths: 15 20 45 20
  :header-rows: 1

  * - Type of tools 
    - Command line
    - Description
    - Default
  * - Pine OMP Threads
    - --pinethreads true
    - Controls thread affinity settings within the OpenMP runtime, influencing how threads are bound to CPU cores for improved performance, particularly on NUMA architectures.
    - false
  * - Set the number of threads
    - --omp_num_threads 10
    - Specifies the number of threads to be utilized for parallel execution, allowing users to control the degree of parallelism based on system resources and workload characteristics.
    - By default it takes the maximum number of threads available
  * - Maximum level of nested parallelism
    - --omp_max_nesting [max_nesting_level]
    - Specifies the maximum level of nested parallelism allowed within OpenMP, controlling the depth at which parallel regions can be nested.
    - -1
  * - Nested parallelism within OpenMP
    - --omp_nested [true/false]
    - Enables or disables nested parallelism within OpenMP, allowing parallel regions to spawn additional parallel regions.
    - false

Tune GPU execution options
--------------------------

Harnessing the power of GPU parallelization, ExaNBody provides users with the ability to fine-tune their execution environment for optimal performance on GPU accelerators. Through a selection of command-line options, users can customize GPU configuration management.

.. list-table:: ExaNBody GPU Command Lines 
  :widths: 15 20 45 20
  :header-rows: 1

  * - Type of tools 
    - Command line
    - Description
    - Default
  * - disable GPU
    - --nogpu
    - disbales use of GPU accelerators, even though some are available.
    - false
  * - workgroup / block size
    - --onika-gpu_block_size N
    - sets default thread block size to N.
    - 128

Profiling tools available in exaNBody
-------------------------------------

ExaNBody offers a comprehensive suite of performance profiling tools designed to empower users in analyzing and optimizing their parallel applications. These tools provide valuable insights into runtime behavior, resource utilization, and performance bottlenecks, enabling developers to fine-tune their code for maximum efficiency. From CPU profiling to memory analysis, ExaNBody's profiling tools offer a range of capabilities to meet diverse profiling needs. This section introduces the profiling tools available within ExaNBody, equipping users with the means to gain deeper understanding and enhance the performance of their parallel applications.

.. list-table:: ExaNBody Profiling Tools Command Lines
  :widths: 15 20 20 45
  :header-rows: 1

  * - Type of tools 
    - Command line
    - Operator
    - Description
  * - Timers 
    - --profiling-summary true
    - profiling : {summary: true }
    - This tool Displays timer informtaions for every operators.
  * - VITE Trace
    - --profilingtrace-file true 
    - NULL
    - This tool generates a VITE trace on CPU (not available with GPU).
  * - Memory footprint 
    - TODO
    - NULL
    - This tool displays the memory footprint of every data storage used during the execution.
  * - nvtx instructions 
    - By default
    - By default
    - Instructions nvtxtoolpush and nvtxtoolpop are included around every operator->execute()
  * - Performance adviser
    - TODO 
    - performance_adviser: { verbose: true }
    - This tool displays some tips according to your simulation (fit cell size, your number of MPI processes ...)


Using Timers with MPI and GPU
------------------------------

In ExaNBody, timers are essential tools for measuring performance in MPI and GPU-accelerated computations. This section explores their use within ExaNBody's parallel implementations, providing insights into runtime behavior and performance characteristics.

This tools provides the list of timers for every operators in a hierarchical form. 
	* Number of calls
	* CPU Time
	* GPU Time
	* Imbalance time between mpi processes (average and maximum)
	* execution time ratio

The Imbalance value is computed as : 
```
I = (T_max - T_ave)/T_ave - 1 
```

With the variables:
	* `T_max` is the execution time of the slowest MPI process.
	* `T_ave` is the average time spent over MPI processes.
	* `I` is the imbalance value.

Note that if you force to stop your simulation, the timer are automatically printed in your terminal.

Output with OpenMP: 

.. code-block:: bash

	Profiling .........................................  tot. time  ( GPU )   avginb  maxinb     count  percent
	sim ...............................................  2.967e+04            0.000   0.000         1  100.00%
	... // some logs
	  loop ............................................  2.964e+04            0.000   0.000         1  99.88%
	    scheme ........................................  2.881e+04            0.000   0.000    100000  97.09%
	      combined_compute_prolog .....................  2.300e+03            0.000   0.000    100000   7.75%
	      check_and_update_particles ..................  1.016e+04            0.000   0.000    100000  34.25%
	        particle_displ_over .......................  2.154e+03            0.000   0.000    100000   7.26%
	        update_particles_full .....................  6.482e+03            0.000   0.000      5961  21.84%
	          update_particles_full_body ..............  6.474e+03            0.000   0.000      5961  21.82%
	            compact_neighbor_friction .............  1.621e+02            0.000   0.000      5961   0.55%
	            move_particles_friction ...............  6.347e+02            0.000   0.000      5961   2.14%
	            trigger_load_balance ..................  2.591e+02            0.000   0.000      5961   0.87%
	              trigger_lb_tmp ......................  6.095e+00            0.000   0.000      5961   0.02%
	                nth_timestep ......................  3.342e+00            0.000   0.000      5961   0.01%
	              extend_domain .......................  2.389e+02            0.000   0.000      5961   0.80%
	...


Output with MPI:

.. code-block:: bash

	Profiling .........................................  tot. time  ( GPU )   avginb  maxinb     count  percent
	sim ...............................................  2.376e+04            0.000   0.000         1  100.00%
	... // some logs
	  loop ............................................  2.372e+04            0.000   0.000         1  99.82%
	    scheme ........................................  2.308e+04            0.086   2.249    100000  97.13%
	      combined_compute_prolog .....................  5.779e+02            0.280   2.937    100000   2.43%
	      check_and_update_particles ..................  1.687e+04            0.454   2.770    100000  70.97%
	        particle_displ_over .......................  4.067e+03            0.687   2.643    100000  17.11%
	        update_particles_full .....................  1.159e+04            0.167   0.812      6001  48.78%
	          update_particles_full_body ..............  1.159e+04            0.167   0.813      6001  48.76%
	            compact_neighbor_friction .............  7.170e+01            0.387   0.876      6001   0.30%
	            move_particles_friction ...............  1.797e+02            0.254   0.853      6001   0.76%
	            trigger_load_balance ..................  9.340e+01            0.674   1.787      6001   0.39%
	              trigger_lb_tmp ......................  2.582e+00            0.187   2.836      6001   0.01%
	                nth_timestep
	              extend_domain .......................  8.655e+01            0.733   2.016      6001   0.36%
	...


Debug features in exaNBody
--------------------------

ExaNBody is equipped with a range of debug features tailored to aid developers in the debugging process. This section outlines the comprehensive list of debug functionalities available within ExaNBody, providing developers with essential tools to diagnose and resolve issues effectively. This is an exhaustive list:

.. list-table:: ExaNBody Debug Command Lines
  :widths: 15 20 20 45
  :header-rows: 1

  * - Type of tools 
    - Command line
    - Architecture
    - Description
  * - Cuda threads size
    - TODO
    - GPU
    - Set the number of cuda threads to 1 on GPU.
  * - Output ldbg
    - --logging-debug true
    - CPU
    - Print debug logs added in `ldbg <<`
  * - filtering debug output
    - --debug-filter ["regexp1","regexp2",...]
    - CPU
    - Filters which operator nodes output debug messges with ldbg<<"...". regexp is a regular expression matching operator pathname, i.e. it's name within block and sub block, for instance "sim.first_iteration.compute_force.lj_force" can be filtered differently than sim.compute_loop.compute_force.lj_force". alternatively, adding a filter expression such as ".*lj_force" will activate debug messages for all instances of lj_force operator.

How to use output ldbg:


Possiblity to active it only for one operator: 
	* Command line : `--logging-debug true --debug-filter[".*operator1",".*operator2",...]`
	* Operator name : logging and debug

Example in your input file (.msp):

.. code-block:: yaml

	configuration:
	  logging: { debug: false , parallel: true }
	  debug:
	    filter: [ ".*init_neighbor_friction" , ".*move_particles_friction" , ".*check_nbh_friction" , ".*compact_neighbor_friction" , ".*extend_domain" ]
