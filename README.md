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
To use the package you can import it using: ```import condor```. Create a ```ccondor_object``` then apply the different functions with ```condor_object.function()```.


## Usage instructions
Supose you have a network (weighted or not) as an edgelist..
```
import condor
co = condor.condor_object("edgelist.txt")
```
Returns a condor object (dictionary with a graph, node names, community membership dataframes...)
```
co.initial_community()
```
Computes the initial community structure and updates the condor object.
```
co.brim()
```
Runs the iterative modularity optimization algorithm (BRIM) and updates the condor object with the final membership.
To see the results type:
```
co.tar_memb
co.reg_memb
```
The whole condor method can be applied to a file by using the function ```run_condor```in the following manner:
```
condor.run_condor("edgelist.txt")
```
it will output ```reg_memb.txt``` and ```tar_memb.txt``` files.

For more detailed information check ```help(condor.condor_object)``` and ```help(condor.run_condor)```.



