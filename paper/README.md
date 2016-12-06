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
    ./prepare_analyses.sh -t $threads -o <mode> [-r $ATROPOS_ROOT]
done

Where <mode> is either 'local' or 'cluster'. This will generate a several command scripts, and a main script, named run_t{threads}_{mode}, that will execute all the necessary steps when run.

On the cluster, we use the 'swarm' script written by Peter Chines (NHGRI), which is specific to Sun Grid Engine. Similar tools should be available for other queue management systems.
