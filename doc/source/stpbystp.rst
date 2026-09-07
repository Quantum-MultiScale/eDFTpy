.. _step_by_step:

===============================
Step by Step Installation Guide
===============================

.. toctree::
   :hidden:
   :maxdepth: 1

   install_types/ubuntu
   install_types/centos
   install_types/amarel


.. raw:: html

    <h2>Before Starting</h2> 


#. **General recommendations:** If you are going to perform the configuration from scratch, we advise carrying out a clean installation on your system. This means deactivating any virtual environment (e.g., Anaconda, Miniconda) before performing the installation to avoid any incompatibility between libraries.

#. **Update Ubuntu/Debian:**
   ::

       $ sudo apt-get update && sudo apt-get upgrade

#. **Update CentOS/RedHat:**
   ::

       $ sudo yum check-update
       $ sudo yum update -y

#. **Check Python version:** Most Unix distributions already include python packages. Before starting your installation, confirm your python version::

       $ python3 --version

#. **Notes about Python >3.10:** If you have Python version > 3.9 (e.g., 3.10.4 or later), you should install libxc from source due to the incompatibility with the pip environment in which eDFTpy was created.

#. **Notes about QE/QEpy compatibility:** eDFTpy and QEpy packages have specific Quantum Espresso dependencies (e.g., version 7.2). Ensure you use a compatible gfortran/gcc compiler for your specific QE build.

#. **Define INSTALL_DIR:** Determine a base directory where you will clone repositories and build dependencies. We will refer to this as your working directory in the following steps.

.. rubric:: Environment Specific Guides

* :doc:`install_types/ubuntu`
* :doc:`install_types/centos`
* :doc:`install_types/amarel`