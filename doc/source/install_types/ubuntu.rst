
.. _ubuntu:

Ubuntu or Debian based OS
=========================

Prerequisites
-------------
Install the necessary system dependencies using the package manager.

#. **cmake, make, git**::

       $ sudo apt-get install cmake make git-all

#. **python3-dev, python3-venv, pip**::

       $ sudo apt-get install python3-dev python3-venv python3-pip

#. **FFTW**::

       $ sudo apt install libfftw3-3 libfftw3-dev

#. **OpenMPI**::

       $ sudo apt-get install libopenmpi-dev openmpi-bin

#. **BLAS/LAPACK**::

       $ sudo apt-get install libblas-dev liblapack-dev

#. **gcc/gfortran**::

       $ sudo apt-get install gcc gfortran

.. note::
   **MKL (optional note):** For maximum performance, you can use Intel MKL libraries (`sudo apt install intel-mkl`) in place of standard BLAS/LAPACK.

Install Packages
----------------

1. **Create virtual environment**::

       $ python3 -m venv edftpy_env
       $ source edftpy_env/bin/activate

2. **Install Python packages**::

       $ pip install --upgrade pip
       $ pip install numpy mpi4py mpi4py-fft wheel setuptools pytest

3. **Compile QE 7.2**::

       $ git clone -b qe-7.2 --depth=1 https://gitlab.com/QEF/q-e.git
       $ cd q-e
       $ ./configure MPIF90=mpif90 --enable-parallel=yes --enable-openmp=no --with-scalapack=no --enable-shared CFLAGS="-fPIC" FFLAGS="-fPIC"
       $ make -j4 all
       $ cd ..
       $ export qedir=$PWD/q-e
       $ export PATH="$qedir/bin:$PATH"
       $ export LD_LIBRARY_PATH="$qedir/lib:$LD_LIBRARY_PATH"

4. **Install QEpy (dev)**::

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