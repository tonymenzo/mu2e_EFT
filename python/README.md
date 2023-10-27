# Introduction

This python code is intended to be a replacment for the Mathematica notebook in ../mathematica
The function names match between the Mathematica and Python code.

The elements data base can be found in elements.yaml in the same directory.

One-body density matrices for relevant  isotopes can be found in a companion repository
[github:Berkeley-Electroweak-Physics/Elastic](https://github.com/Berkeley-Electroweak-Physics/Elastic) .  You should clone this repository and set the environment variable MU2E_ELASTIC to it's full path.
If you were added to the access list for the Mu2E scripts then you should be on the access list for the Elastic directory.   If not, just ask.
```
export MU2E_ELASTIC=<clonepath>/Elastic
```

# Running
Your input setup should be in a <mycase>.yaml file.

```
python3 mu2e.py [args]
python3 mu2e.py -help   # generates help message
```

A standard run will consume a yaml file specifying the required
data for the run.   "example.yaml" in this directory is a heavily commented template for 63Cu.  
```
python3 mu2e.py example.yaml  # run file example.yaml
```

mu2e also has an option to produce a Mathematica script to run the Mathematica version of the code.
```
python3 mu2e.py -wl example.yaml
```
This run will additionally produce example.wl (wl is a standard suffix for Mathematica files).    The beginning of the file will
explain how to run it.


The examples subdirectory has additional isotope specific examples.

Plots will be written to a directory who's name is the name of the element being processed.

# Requirements

```
pip3 install argparse         # command line parsing
pip3 install pathlib          # File system api
pip3 install multipledispatch # function overloading
pip3 install pyyaml           # python data files.  
                              # Used for isotope data
pip3 install py3nj            # implements 3J, 6J, 9J symbols
pip3 install numpy            # arrays, matrices, tensors, ...
pip3 install scipy            # A variety of special functions
                              #  hypergeometric, integration, ...
pip3 install matplotlib       # Powerful plotting package.
```
