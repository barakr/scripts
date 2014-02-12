#!/bin/bash
#$ -S /bin/bash
#$ -o /scrapp/barak/Logs/
#$ -cwd
#$ -j y
#$ -r y
#$ -N cyl.job
#$ -l arch=linux-x64,mem_free=1.5G
#$ -l h_rt=1:00:00
#$ -t 1

#if [ $# -lt 2 ] ; then
#    echo "Usage: $0 <flags_file> <serial number>"
#    exit -1
#fi
echo Running \"$0 $*\"
export DECOYS=decoys_flags_abinitio.1.silent #decoys_$1.$2.silent
if [ ! -e $DECOYS ]; then
    echo "file $DECOYS does not exist";
    exit $STATUS;
fi
./run_src/Scripts/cluster.sh 500 2 score_flags_abinitio.1.sc run_src/native.pdb $DECOYS reweighted_sc
