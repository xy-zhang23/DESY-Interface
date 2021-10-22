# README #


### What is this repository for? ###

* An interface script for **Astra** and **Genesis**
* Version 1.0.0

### How do I get set up? ###

Download the repository

```bash
git clone https://XiangkunLi@bitbucket.org/XiangkunLi/interface.git
```

Go to the folder where `setup.py` is and run

```bash
python setup.py install
```

To load the modules, use

```python
from interface import *
```


### Tutorials with Jupyter notebook ###
- For `Astra` related interface, see `tutorials/astra_demo.ipynb`  (or `tutorials/astra_demo.html` in `html` format)
- For batch generating Genesis 1.3 input files, see `tutorials/genesis13_demo.ipynb'
- For postprocessing of Genesis 1.3 simulations, see `tutorials/postG4_demo.ipynb`


### Examples ###

Examples include:

- `dogleg_demo.py`, set up a dogleg for astra simulation
- `injector_demo.py`, set up a photo-injector for astra simulation
- `injector_optimization_demo.py` and `sub.sh`, set up an optimization script using MOGA algorithm and submit the job to cluster.

About the `NSGA` algorithm , the default is the python module `platypus`. The problem with it is that there is no intermediate 
output and one simply doesn't know what happened in the optimization procedure. Therefore, a wrapper class `interface/NSGAPlus` was 
written to solve this problem. The class is just modified from the default one, with some options to control the output of the intermediate results.

### About the `Namelist/Generator1/Astra` class ###

There are mainly three data types involved to use these classes:

- non-array type, e.g., int, float, double or boolean: used to define variables, such as `sig_x` or `Zstop1`
- 1D array-like type, e.g., a tuple: (1, 2, 3), or a list [1, 2, 3], or an `Numpy` array: numpy.array([1, 2, 3]), used to define 1D variables, such as `MaxE`
- 2D array-like type, could also be tuple or list or `Numpy` array, used to define special 2D variables such as `D1` from the `Dipole` namelist

The fieldmaps are define with 1D array-like type, the elements of which are strings (`str` type). 
Note that the quadrupoles and 3D fieldmaps seem not support the use of absolute paths.

### Bugs and updates

*02.10.2021*

- Renamed **Namelist.py** to **Namelists.p**; renamed **postG4.py** to **postGenesis13.py**. The relevant importing in other files should be updated as well.