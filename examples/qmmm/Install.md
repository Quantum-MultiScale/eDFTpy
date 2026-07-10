# Installation in Amareln

#### 1. Load modules

```shell
module use /projects/community/modulefiles
module load python/3.8.5-gc563 intel/19.1.1
```

#### 2. Create an environment 

```shell
python3 -m venv edftpy_qmmm in software folder
source ~/edftpy_qmmm/bin/activate
```

#### 3. Export some libraries

```shell
- export LD_PRELOAD=/opt/sw/packages/intel/19.1.3/mkl/lib/intel64/libmkl_rt.so
- export LD_LIBRARY_PATH="/opt/sw/packages/intel/19.1.3/mkl/lib/intel64/:$LD_LIBRARY_PATH"
- export LDFLAGS="-L//opt/sw/packages/intel/19.1.3/mkl/lib/intel64"
```

#### 4. Installation of QE

```shell
git clone https://gitlab.com/QEF/q-e.git
git checkout qe-7.2
git submodule update --recursive
module load make/4.2-kholodvl git
make distclean
export IFLAGS="-I~/q-e/external/devxlib/src -I~/q-e/include -I~/q-e/external/devxlib/include/"
- ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC 
- make all -j6
- export qedir=`pwd`
```

#### 4. Install qepy 

```shell
pip install ninja
pip install f90wrap
pip install mpi4py-fft pytest pytest-cov pylibxc2 pytest-mpi
git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
cd qepy
git checkout dev
tddft=yes python -m pip install -U .
```

#### 5. Install DFTpy 

```shell
python -m pip install scipy
pip install ase==3.22.0
pip install xmltodict
pip install -i https://test.pypi.org/simple/ dftpy==2.1.0 
```

OR

```shell
git clone https://gitlab.com/pavanello-research-group/dftpy.git
python -m pip install .
```

#### 6. Install eDFTpy version qmmm

```shell
pip install wheel
python -m pip install --upgrade setuptools wheel
git clone https://gitlab.com/pavanello-research-group/edftpy.git
git checkout qmmm-dev
python -m pip install .
```

#### 7.  Install MBX
```shell
module load gsl/2.5-bd387
git clone https://github.com/paesanilab/MBX.git (Do not use it)
git clone https://gitlab.com/xinchen_compchem/mbx-python-old.git
cd mbx-python-old
autoreconf -fi
./configure --disable-optimization --enable-shared --prefix=~/software/mbx-python-old/install/
```

##### 7.1 Modify make file

- Change : LIBS = -lgslcblas -lgsl -lfftw3  --> To LIBS = -lgsl -lgslcblas -lfftw3 -liomp5  -lpthread -fopenmp 

##### 7.2 Compile
```shell
- make -j6 && make install 
```

####  8. Install Kernel
```shell
pip install ipykernel
ipython kernel install --user --name=edftpy_qmmm
```

#### 9. Add MBX to python path

```shell
export PYTHONPATH=`$`PYTHONPATH:`$`{MBX_HOME}/plugins/python/mbx
```

#### 10. Install DFT-D3
```shell
module load gcc/
meson setup _build --prefix=~/software/
meson test -C _build --print-errorlogs
meson install -C _build

module load gcc/10.2.0/openmpi
export LD_LIBRARY_PATH="~/software/dftd3/lib64/:$LD_LIBRARY_PATH"
dftd3='~/software/dftd3/bin/s-dftd3'
```

