#!/bin/bash
#$ -S /bin/bash
#$ -l arch=linux-x64,mem_free=0.75G
#$ -l h_rt=72:00:00
#$ -pe ompi 100
#$ -R yes
#$ -V
#$ -cwd
#$ -o /scrapp/barak/Logs
#$ -N FPDock

module load openmpi-x86_64

if [ $# -lt 3 ] ; then
    echo "Usage: $0 <flags_file> <serial number> <out_folder_suffix>"
    exit -1
fi
echo Running \"$0 $*\"

# Prepare
SRC=`pwd -P`
echo $SRC
echo WHO
OUTFOLDER=/scrapp/barak/FPDock_$3
echo $OUTFOLDER
if [ ! -d $OUTFOLDER ]; then
    mkdir $OUTFOLDER;
    STATUS=$?
    if [ $STATUS != 0 ]; then echo "couldn't create $OUTFOLDER"; exit $STATUS; fi
fi
cd $OUTFOLDER
echo HO
pwd -P
ln -s $SRC Run
#set MYTMP=`mktemp -d`
#set seed=`od -An -N4 -td4 /dev/random`

# Run:
#cd $MYTMP
#echo "Temporary run folder $MYTMP"
mpirun -np $NSLOTS /netapp/sali/barak/Rosetta/rosetta_source/bin/FlexPepDocking.mpi.linuxgccrelease -database /netapp/sali/barak/Rosetta/rosetta_database/  @Run/$1 -out:file:silent decoys_$1.$2.silent -scorefile score_$1.$2.sc

echo HI

Run/Scripts/SubmitAndCluster/cluster.sh 500 2 score_$1.$2.sc Run/native.pdb decoys_$1.$2.silent reweighted_sc

#echo "Moving final output files from $MYTMP to $OUTFOLDER"
#mv $MYTMP/* $OUTFOLDER
#rmdir $MYTMP
