The scripts in this directory should enable you to re-run the analyses in the Atropos paper.

First, look in the header of the install_software.sh script to make sure you have the prerequisites installed. Note that atropos will be installed via pip, and then added to the current virtual environment. You can create and use a new virtual environment for testing purposes - and this is required if your normal environment has python < 3.3:

  conda create --name atropos_test python=3

Next, install software needed to simulate reads, and all of the tools we will benchmark:

  export ATROPOS_ROOT=/path/to/atropos/paper
  cd $ATROPOS_ROOT/scripts
  ./install_software.sh

Simulate reads:
  
  ./simulate_reads.sh

Prepare the command scripts:

  for threads in 4 8 16
  do
    ./prepare_analyses.sh -t $threads -r $ATROPOS_ROOT
  done

We executed the commands both locally (on a desktop computer) and on a remote cluster. On the cluster, we use the 'swarm' script written by Peter Chines (NHGRI), which is specific to Sun Grid Engine. Similar tools should be available for other queue management systems.
  
  # on desktop
  for threads in 4 8 16
  do
    ./commands_t${threads}.sh 2> timing_${threads}
    # summarize timing info
    python summarize_timing_info.py -i timing_${threads} -o ../results/local_timing_t${threads}.txt
  done
  
  # on cluster
  for threads in 4 8 16
  do
    swarm --threads-per-process 16 --gb-per-process 4 --file commands_t${threads}.sh --name atropos_test_t${threads}
    # SGE generates a .eXXX and .oXXX file for each process. Timing info is in the .eXXX file.
    qsub -hold_jid atropos_test_t${threads} -b y \
      cat commands_t${threads}.e* | python summarize_timing_info.py -o ../results/cluster_timing_t${threads}.txt
  done

Finally, we compute the accuracy benchmarks. Since (theoretically) the results should be exactly the same regardless of the number of threads used, we arbitrarily choose the 4 threads results to analyze.

  for err in 001 005 01
  do
    for qcut in 0 20
    do
      python summarize_accuracy.py
