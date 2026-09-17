.. _amarel:

Amarel Cluster Installation
===========================

.. note::
    Installation of eDFTpy on Amarel from Scratch (In a place where we can all access it!)


Prerequisites
-------------

.. raw:: html

    <h2>Required modules</h2>

Load the necessary modules for the compilation environment::

    $ module use /projects/community-old/modulefiles
    $ module load gcc/15.2.0/openssl/3.6.0-krasting libffi/3.3-gc563
    $ module load intel/19.1.1
    $ module load gcc/14.2.0-cermak

.. raw:: html

    <h2>Define INSTALL_DIR</h2>

Set your installation prefixes::

    $ export PREFIX=/home/user/where_to_install/python310
    $ export PREFIX2=/home/user/where_to_install/

Install Packages
----------------

.. raw:: html

    <h2>Compile Python 3.10</h2>

Compile Python to ensure full compatibility with the loaded modules::

    $ export CPPFLAGS="-I/projects/community-old/libffi/3.3/gc563/include -I/projects/community-old/gcc/15.2.0/openssl/3.6.0-krasting/include"
    $ export LDFLAGS="-L/projects/community-old/libffi/3.3/gc563/lib64 -L/projects/community-old/gcc/15.2.0/openssl/3.6.0-krasting/lib64"
    $ export LIBS="-lffi -lssl -lcrypto"

    $ mkdir -p $PREFIX
    $ ./configure --prefix=$PREFIX --enable-optimizations --with-ensurepip=install
    $ make -j4
    $ make install

Verify the Python installation::

    $ python3.10 -c "import _ctypes; import ssl; print('OK')"
    $ export PATH=$PREFIX/bin:$PATH

.. raw:: html

    <h2>Create virtual environment</h2>

Create and activate the environment::

    $ python3.10 -m venv edftpy_2026_2
    $ source "$PREFIX2/edftpy_2026_2/bin/activate"
    $ pip install --upgrade pip

.. raw:: html

    <h2>Compile FFTW</h2>

Build FFTW with Intel compiler support::

    $ wget http://www.fftw.org/fftw-3.3.10.tar.gz
    $ tar -xzf fftw-3.3.10.tar.gz
    $ cd fftw-3.3.10
    $ ./configure --prefix=/home/user/where_to_install/fftw_intel --enable-shared --enable-threads --with-pic --enable-mpi
    $ make -j4
    $ make install
    $ cd ../fftw_intel

    $ export FFTW_HOME=/home/user/where_to_install/fftw_intel
    $ export FFTW_DIR=$FFTW_HOME
    $ export FFTW_LIB_DIR=$FFTW_DIR/lib
    $ export FFTW_INCLUDE_DIR=$FFTW_DIR/include
    $ export LD_LIBRARY_PATH=$FFTW_HOME/lib:$LD_LIBRARY_PATH

.. raw:: html

    <h2>Install mpi4py/mpi4py-fft</h2>

Export compiler flags and install the required Python wrappers::

    $ export CC=mpicc
    $ export CXX=mpicxx
    $ export MPICC=mpicc
    $ export CFLAGS="-std=gnu99"
    $ export CXXFLAGS="-std=c++11"

    $ pip install mpi4py
    $ pip install mpi4py-fft
    $ python -m pip install wheel pytest pytest-cov pylibxc2 pytest-mpi

.. raw:: html

    <h2>Compile QE 7.2 (MKL)</h2>

Install a parallel version of QEpy that is compatible with all features of eDFTpy::

    $ git clone -b qe-7.2 --depth=1 https://gitlab.com/QEF/q-e.git
    $ cd q-e/
    $ git submodule update --recursive
    $ python -m pip install --upgrade setuptools wheel

    $ export IFLAGS="-I$PWD/include -I$PWD/external/devxlib/src -I$PWD/external/devxlib/include"
    $ export BLAS_LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl"
    $ export LAPACK_LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential"

    $ ./configure MPIF90=$MPIF90 --enable-parallel=yes --enable-openmp=no --with-scalapack=no --enable-shared CFLAGS="-fPIC" FFLAGS="-fPIC"
    $ make -j4 all

    $ QEDIR=/home/user/where_to_install/q-e
    $ export qedir=$QEDIR
    $ export PATH="$QEDIR/bin:$PATH"
    $ export LD_LIBRARY_PATH="$QEDIR/lib:$LD_LIBRARY_PATH"

.. raw:: html

    <h2>Install QEpy (dev)</h2>

::

    $ git clone --recurse-submodules https://github.com/Quantum-MultiScale/QEpy.git
    $ cd QEpy
    $ git checkout dev
    $ tddft=yes python -m pip install -U .
    $ cd ..

.. raw:: html

    <h2>Install DFTpy</h2>

::

    $ git clone --recurse-submodules https://github.com/Quantum-MultiScale/DFTpy.git
    $ cd DFTpy
    $ # (Not building from source for now until Michele's bug is solved)
    $ pip install dftpy==2.2.0
    $ cd ..

.. raw:: html

    <h2>Install eDFTpy</h2>

::

    $ git clone --recurse-submodules https://github.com/Quantum-MultiScale/eDFTpy.git
    $ cd eDFTpy
    $ git checkout dev
    $ python -m pip install .
    $ cd ..

.. raw:: html

    <h2>Compile LibXC</h2>

::

    $ module load cmake/3.24.3-sw1088
    $ git clone https://gitlab.com/libxc/libxc.git
    $ cd libxc
    $ export LIBRARY_PATH=/projects/community-old/lapack/3.11.0/sw1088/lib64:$LIBRARY_PATH
    $ export LD_LIBRARY_PATH=/projects/community-old/lapack/3.11.0/sw1088/lib64:$LD_LIBRARY_PATH
    $ mkdir build
    $ cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DDISABLE_KXC=OFF -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$PWD/build/lib
    $ cmake --build build -j6
    $ cmake --install build --prefix=/home/user/where_to_install/libxc
    $ cd ..

.. raw:: html

    <h2>Verify installation</h2>

::

    $ python -c "import edftpy; import qepy; print('Installation Successful!')"