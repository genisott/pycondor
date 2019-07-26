# pyCONDOR
Python implementation of the BRIM algorithm for bipartite community structure detection as described in "Modularity and community detection in bipartite networks" by Michael J. Barber." This package is somewhat translated from the R package CONDOR https://github.com/jplatig/condor described in the paper "Bipartite Community Structure of eQTLs" by John Platig , Peter J. Castaldi, Dawn DeMeo, John Quackenbush.

## Install instructions
This package uses the python libraries numpy,pandas and python-igraph (igraph).

The easiest way to install pyCONDOR is using pip:
```
git clone https://github.com/genisott/pycondor.git
cd pycondor
pip install .
```
To use the package you can import it using: ```import condor```. And then access the different functions implemented with ```condor.function()```.
