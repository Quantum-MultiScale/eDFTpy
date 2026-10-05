#bash

source /home/jezs/Documents/4_Embedding/edftpy_2026/bin/activate

export LD_LIBRARY_PATH=/home/jezs/Documents/3_QMMM/software/libxc/build/lib.linux-x86_64-cpython-310/pylibxc:$LD_LIBRARY_PATH

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "Python Path:"

python -c "import pylibxc; print('pylibxc OK — core path:', getattr(pylibxc,'get_core_path',lambda:None)())"

#export MBX_HOME=/home/jezs/Documents/3_QMMM/software/mbx-python-old/
export MBX_HOME=/home/jezs/Documents/3_QMMM/new_software/mbx-dipoles/
export PYTHONPATH=$PYTHONPATH:${MBX_HOME}/plugins/python/mbx

mpirun -n 4 python test_QMMM.py > log
