#bash

source /home/jezs/Documents/4_Embedding/edftpy_2026/bin/activate

export LD_LIBRARY_PATH=/home/jezs/Documents/3_QMMM/software/libxc/build/lib.linux-x86_64-cpython-310/pylibxc:$LD_LIBRARY_PATH

if [ -z "$MPIRUN" ]
then
	mpirun="mpirun -n 4"
else
	mpirun=$MPIRUN
fi
	
$mpirun python -m edftpy edftpy_1.ini | tee log.1
