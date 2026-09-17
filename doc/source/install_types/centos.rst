
.. _centos:

CentOS and RedHat based OS
==========================

Prerequisites
-------------

#. **yum packages and python3**::

       $ sudo yum install epel-release
       $ sudo yum update -y
       $ sudo yum install -y python3 python3-devel python3-pip

#. **FFTW**::

       $ sudo yum install -y fftw-devel

#. **OpenMPI**::

       $ sudo yum install -y openmpi openmpi-devel
       $ module load mpi/openmpi-x86_64

#. **gcc/gfortran**::

       $ sudo yum install -y gcc gcc-gfortran gcc-c++

#. **BLAS/LAPACK**::

       $ sudo yum install -y lapack-devel blas-devel

#. **cmake and git**::

       $ sudo yum install -y cmake make git

Install Packages
----------------

1. **Create virtual environment**::

       $ python3 -m venv edftpy_env
       $ source edftpy_env/bin/activate

2. **Install Python packages**::

       $ pip install --upgrade pip
       $ pip install numpy cython mpi4py mpi4py-fft wheel setuptools

3. **Compile QE**::

       $ git clone -b qe-7.2 --depth=1 https://gitlab.com/QEF/q-e.git
       $ cd q-e
       $ ./configure MPIF90=mpif90 --enable-parallel=yes --enable-openmp=no --with-scalapack=no --enable-shared CFLAGS="-fPIC" FFLAGS="-fPIC"
       $ make -j4 all
       $ cd ..
       $ export qedir=$PWD/q-e
       $ export PATH="$qedir/bin:$PATH"
       $ export LD_LIBRARY_PATH="$qedir/lib:$LD_LIBRARY_PATH"

4. **Install QEpy**::

       $ git clone --recurse-submodules https://github.com/Quantum-MultiScale/QEpy.git
       $ cd QEpy
       $ git checkout dev
       $ tddft=yes python -m pip install -U .
       $ cd ..

5. **Install DFTpy**::

       $ pip install dftpy==2.2.0

6. **Install eDFTpy**::

       $ git clone --recurse-submodules https://github.com/Quantum-MultiScale/eDFTpy.git
       $ cd eDFTpy
       $ git checkout dev
       $ python -m pip install .
       $ cd ..

7. **Compile LibXC**::

       $ git clone https://gitlab.com/libxc/libxc.git
       $ cd libxc
       $ mkdir build
       $ cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DDISABLE_KXC=OFF -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$PWD/build/lib
       $ cmake --build build -j4
       $ sudo cmake --install build
       $ cd ..

8. **Verify installation**::

       $ python -c "import edftpy; import qepy; print('Installation Successful!')"