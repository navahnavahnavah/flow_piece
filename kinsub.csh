#!/bin/tcsh
#
# richardd, 7 Oct 11
#
# converted to tcsh, 17 Oct 11
#
# example PBS script for berserker
#
# put this script in the same directory as the program
# you want to run.
#
# set the name of the job
#PBS -N fp_550_250_pres_q1
#
# set the output and error files
#PBS -o /data/navah/fp_output/$PBS_JOBNAME/mOutG.txt
#PBS -e /data/navah/fp_output/$PBS_JOBNAME/mErrG.txt
#PBS -m abe -M navah@uchicago.edu
# set the number of nodes to use, and number of processors
# to use per node


#PBS -l nodes=compute-1-3:ppn=1


# or, if using only one node, you can do it this way too
##PBS -l ncpus=5


# NEW STUFF MAY 2016

# in this example, I'm using the intel compilers and mvapich2
#
# bring in the module settings

sleep 10s
source /etc/profile.d/modules.csh
module load intel/intel-12
module load mpi/mvapich2/intel

set PARAM_RESTART='0'
set PARAM_PATH='/data/navah/fp_output/'$PBS_JOBNAME'/'
set PARAM_CRASHSTEP='2'


# model parameters go here i guess
set PARAM_O='550'
set PARAM_W='300'
set PARAM_W_RHS='25'
set PARAM_H='200'
# set PARAM_O_RHS='-50'
set PARAM_O_RHS='0'
set PARAM_H_S='250'






set PARAM_PAQ = '1e-12'

## FRACTURE PARAMS
set PARAM_F_DX = '-1.1'
set PARAM_F_POR = '2.5e-5'








echo $PARAM_PATH
# set PROGNAME to the name of your program
set PROGNAME=massacr

# figure out which mpiexec to use
set LAUNCH=/usr/mpi/intel/openmpi-1.4.3-qlc/bin/mpirun

# working directory
set WORKDIR=${HOME}
set WORKDIR=/home/navah/flow_piece

set NCPU=`wc -l < $PBS_NODEFILE`
set NNODES=`uniq $PBS_NODEFILE | wc -l`
# set this to zero to turn OFF debugging, 1 to turn it on
set PERMDIR=${HOME}
set SERVPERMDIR=${PBS_O_HOST}:${PERMDIR}

set DEBUG=1
if ( $DEBUG ) then
	echo ------------------------------------------------------
	echo ' This job is allocated on '${NCPU}' cpu(s)'
	echo 'Job is running on node(s): '
	cat $PBS_NODEFILE
	echo ------------------------------------------------------
	echo PBS: qsub is running on $PBS_O_HOST
	echo PBS: originating queue is $PBS_O_QUEUE
	echo PBS: executing queue is $PBS_QUEUE
	echo PBS: working directory is $PBS_O_WORKDIR
	echo PBS: execution mode is $PBS_ENVIRONMENT
	echo PBS: job identifier is $PBS_JOBID
	echo PBS: job name is $PBS_JOBNAME
	echo PBS: node file is $PBS_NODEFILE
	echo PBS: number of nodes is $NNODES
	echo PBS: current home directory is $PBS_O_HOME
	echo PBS: PATH = $PBS_O_PATH
	echo ------------------------------------------------------
	echo workdir is $WORKDIR
	echo permdir is $PERMDIR
	echo servpermdir is $SERVPERMDIR
	echo ------------------------------------------------------
	echo 'Job is running on node(s): '
	cat $PBS_NODEFILE
	echo ------------------------------------------------------
	echo ${LAUNCH} -n {$NCPU} -f ${PBS_NODEFILE} ${WORKDIR}/${PROGNAME}
	echo ' '
	echo ' '
endif

echo $NCPU
echo $PBS_NODEFILE
echo $PARAM_PATH

# run the program

# set SCRDIR=/state/partition1/navah/${PBS_JOBID}/
# mkdir -p $SCRDIR


#cd $SCRDIR
cd ${WORKDIR}
##wait
${LAUNCH} -n {$NCPU} -hostfile ${PBS_NODEFILE} --nooversubscribe ${WORKDIR}/${PROGNAME} ${PARAM_RESTART} ${PARAM_PATH} ${PARAM_CRASHSTEP} ${PARAM_O} ${PARAM_W} ${PARAM_W_RHS} ${PARAM_H} ${PARAM_O_RHS} ${PARAM_PAQ} ${PARAM_F_DX} ${PARAM_F_POR} ${PARAM_H_S}


# wait
# ssh $PBS_NODEFILE
# cd $SCRDIR
# mv * $PARAM_PATH_ALT
