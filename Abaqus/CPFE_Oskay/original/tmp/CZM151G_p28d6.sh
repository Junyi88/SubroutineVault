# Job Name
#$ -N CZM151G_p28d6
# Use current working directory
#$ -cwd

# Parallel Environment request.  Set your number of processors here
#$ -pe mpich 16

# Run job through bash shell
#$ -S /bin/bash


# Add the intel Fortran Compiler Module

/cm/shared/apps/intel/parallel_studio_xe/current/bin/ifortvars.sh intel64

# Tell SGE to show the PATH variable

PATH=$PATH:/cm/shared/licenses/intel:/cm/shared/apps/intel/parallel_studio_xe/current/bin/intel64
PATH=$PATH:/cm/shared/apps/SIMULIA/6.14-3:/cm/shared/apps/SIMULIA/6.14-3/code/Python2.7:/cm/shared/apps/SIMULIA/6.14-3/Python2.7/Lib

LD_LIBRARY_PATH=/cm/shared/apps/intel/parallel_studio_xe/current/compiler/lib/intel64


echo Got $NSLOTS processors.
echo Machines:
echo JOBNAME: $JOBNAME
echo JOBN0: $JOBN0
echo JOB_ID: ${JOB_ID}
cat $TMPDIR/machines


# Define the Abaqus particulars
JOB_NAME=CZM151G_p28d6
INPUT_FILE=/home/tphan/tmp/${JOB_NAME}.inp

ABAQUS_ARGS="cpus=16 mp_mode=mpi user=/home/tphan/tmp/umat_uel.f -inter"
ABAQUS="/cm/shared/apps/SIMULIA/Commands/abaqus"

# get the current working directory
cwd=$(pwd)
echo pwd: $cwd
cd $cwd/
## Run the job
rm ${JOB_NAME}.lck
$ABAQUS job=${JOB_NAME} input=${INPUT_FILE} ${ABAQUS_ARGS}
# $ABAQUS python ReadData.py
# rm teture.txto
