# markedgraphs
### A Python module for working with marked graphs.
By Stefan van der Lugt (Universiteit Leiden). 
An implementation of the algorithms in my PhD thesis *Tautological differential forms on moduli spaces of curves*.

## The class `MarkedGraph`
The `markedgraphs` module provides a class `MarkedGraph` that serves as a model for marked graphs.

A marked graph can be constructed in two ways: by providing a list of edges and by providing an adjacency matrix.
For example, creating the 1-marked graph with 1 unmarked vertices and 3 edges between these vertices can be done in the following ways:
```python
from markedgraphs import MarkedGraph
G1 = MarkedGraph(
        marked=1, 
        vertices=2, 
        edges=[(0,1),(0,1),(0,1)])

G2 = MarkedGraph(
    marked = 1,
    mat = [[0,3],[3,0]])
```
When giving an adjacency matrix, we assume that the first rows/columns are the ones corresponding to the marked vertices.

In order to check if two graphs are isomorphic one can check for equality.
For instance, the boolean `G1 == G2` is `True` in the above example.

The class `MarkedGraphs` has various methods for manipulating marked graphs.
An example:
```python
from markedgraphs import MarkedGraph

G3 = MarkedGraph(
    marked=0,
    vertices=2,
    edges=[(0,0),(0,1),(1,1)])

G4 = MarkedGraph(
    marked=0,
    vertices=1,
    edges=[(0,0),(0,0)])

print(G3.is_contracted()) # False
print(G4.is_contracted()) # True
print(G3 == G4)           # False
G3.contract()
print(G3.is_contracted()) # True
print(G3==G4)             # True 
```

## The function `CG`
The module also provides a function `CG(r,chi,u=None)` that returns a list of all contracted *r*-marked graphs of characteristic *χ* with *u* unmarked vertices.
If *u* is not given, `CG(r,chi)` returns a list of all contracted *r*-marked graphs of characteristic *χ*.

Example: the following script prints the adjacency matrices of all 21 contracted 1-marked graphs of characteristic -1:
```python
from markedgraphs import CG

for G in CG(1,-1):
    print(G.matrix(),"\n")
```

## Example scripts
The `examples/` folder contains some scripts that were used in my thesis.

The script `compute_fd.py` allows the user to compute the polynomial *f<sub>d,u</sub>(r)* for all possible values of *d* and *u*.
This polynomial gives a closed expression for the number of contracted *r*-marked graphs of characteristic *r-d* with *u* unmarked vertices.
Likewise, the polynomial *f<sub>d</sub>(r)* that gives the number of contracted *r*-marked graphs of characteristic *r-d* can be computed by not providing any value for *u*. 

The script `find_relations.py` can be used to compute some relations among tautological differential forms by using the identity *ω<sub>0</sub><sup>g+1</sup>=0*.
It asks for integers r, s, and G, finds some relations in *R<sup>2G+2</sup>(C<sub>g</sub><sup>r</sup>)*, and takes fiber integrals to obtain relations in *R<sup>2G+2-2r+2s</sup>(C<sub>g</sub><sup>s</sup>)*.

## Requirements
The `markedgraphs` module requires `numpy`.

The script `compute_fd.py` needs `sympy` to carry out Lagrange interpolation.
The script `find_relations.py` uses `sagemath` to work with polynomial rings and vector spaces.
