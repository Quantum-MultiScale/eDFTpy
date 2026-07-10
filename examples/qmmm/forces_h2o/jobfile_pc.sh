#bash

installation='new_software'
env_name='edftpy_qmmm_new'

#installation='software'
#env_name='edftpy_qmmm'

env=~/Documents/3_QMMM/$installation/$env_name/bin/activate

mbx='mbx-dipoles'
#mbx='mbx-no-dipole'

export LD_LIBRARY_PATH=/home/jezs/Documents/3_QMMM/software/libxc/build/lib.linux-x86_64-cpython-310/pylibxc:$LD_LIBRARY_PATH

source $env
export MBX_HOME=/home/jezs/Documents/3_QMMM/$installation/$mbx/
export LD_LIBRARY_PATH=/home/jezs/Documents/3_QMMM/$installation/$mbx/install/lib
export PYTHONPATH=$PYTHONPATH:${MBX_HOME}/plugins/python/mbx

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "Python Path:"
python -c 'import sys; print(sys.path)'

for ifile in xxyz06--2.xyz xxyz06--1.xyz xxyz06.xyz xxyz06-+1.xyz xxyz06-+2.xyz ;
do
  echo "$ifile"
  #mv $ifile $ifile'.xyz'
  sed "19c cell-file = $ifile"  mbx-tmp.ini > mbx.ini
#  srun -N $NODES -n $NTASKS --mpi=pmi2 python3 -m edftpy --run mbx.ini --mpi4py  |tee $ifile.out
  mpirun -np 3 python -m edftpy --run mbx.ini --mpi4py |tee $ifile.out
  wait
  cp  sub_ks1.out   $ifile'_ks1.rec' 
  cp  sub_mm.out   $ifile'_mm.rec' 
done

