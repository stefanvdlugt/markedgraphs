# This script is a quicker version of find_relations.py
# Rather than a CombinatorialFreeModule it uses a VectorSpace which is faster when it comes to computing bases.

from sage.all import FractionField, PolynomialRing, QQ, VectorSpace

# allow import of markedgraphs module from parent directory:
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from markedgraphs import CG, MarkedGraph
from collections import defaultdict
from math import gcd
from functools import reduce

# The following function implements a cache dictionary of all graphs we encountered
seen_graphs = defaultdict(list) # Cache of all graphs encountered before
def graph_index(Gamma):
    """
    Takes a contracted marked graph G, and return a tuple (r,chi,i) where:
    - r is the number of marked vertices of G;
    - chi is the characteristic of G;
    - i is a nonnegative integer that indicates the isomorphism class of G.
    Warning: i is not an invariant of G, it may change between runs 
    """
    r = Gamma.marked()
    chi = Gamma.characteristic()
    for i,Gamma2 in enumerate(seen_graphs[r,chi]):
        if Gamma == Gamma2:
            return (r,chi,i)
    else:
        seen_graphs[r,chi].append(Gamma)
        return (r,chi,len(seen_graphs[r,chi])-1)

def lcm(a,b):
    return a*b//gcd(a,b)

r = int(input('Enter r: '))
s = int(input('Enter s: '))
G = int(input('Enter G: '))

d = 2*(G+1-r+s) # degree of resulting forms on C_g^s

# Construct the polynomial ring Q[g]
QQg = PolynomialRing(QQ,'g')
g = QQg.gens_dict()['g']

# Construct the polynomial ring Z[g][Gamma_{i,j}: 1 <= i <= j <= r].
# There is a canonical homomorphism Z[g][Gamma_{i,j}] --> Z[g][CG(r)] 
# that maps Gamma_{i,j} to Gamma_{i,j}.
Rvars = [f'Gamma_{i}_{j}' for i in range(r) for j in range(i,r)]
R = QQg[','.join(Rvars)]
Gamma = {(i,j): R.gens_dict()[f'Gamma_{min(i,j)}_{max(i,j)}'] for i in range(r) for j in range(r)}

# Construct the polynomial ring Z[g][Gamma_{i,j}][x_k] with 1<=k<r.
Svars = [f'x_{i}' for i in range(r-1)]
S = R[','.join(Svars)]
x = {i: S.gens_dict()[f'x_{i}'] for i in range(r-1)}
x[r-1] = -sum(x.values())

# Compute W_r
Wr = sum(Gamma[i,j]*x[i]*x[j] for i in range(r) for j in range(r))

# Compute W_r^{G+1} and take its coefficients; these are elements of Z[g][Gamma_{ij}]
# and contained in the kernel of Z[g][Gamma_{ij}] --> Z[g][CG(r)] --> R(C_g^r)
relations = (Wr**(G+1)).coefficients()

# For each of the coefficients found, take the image 
# under the homomorphism Z[g][Gamma_{i,j}] --> Z[g][CG(r)] --> Z[g][G(s)] --> Z[g][CG(s)]
# We'll represent the resulting elements of Z[g][CG(s)] as dicts. 
# The keys correspond to elements of CG(s), and its values are the coefficients in Z[g].

# edge[k] = Gamma_{i,j} iff Gamma_{i,j} is the k-th generator of R.
edge = {R.gens().index(Gamma[i,j]): (i,j) for i in range(r) for j in range(i,r)}
allkeys = set() # Lists graphs with nonzero coefficient
CGs_elements = list()
for rel in relations:
    elt = defaultdict(QQg)
    for tup, coeff in rel.dict().items():
        Gr = MarkedGraph(marked=r, vertices=r, edges=[])
        for k in range(len(tup)):
            for _ in range(tup[k]):
                Gr.add_edge(*edge[k])
        Gs = MarkedGraph(marked=s, mat = Gr.matrix())
        lam = Gs.contract(g=g)
        Gs_index = graph_index(Gs)
        elt[Gs_index] += lam*coeff
    
    # remove terms with coefficient 0
    elt = {k:v for k,v in elt.items() if v!=0}
    if elt:
        CGs_elements.append(elt)
        allkeys |= elt.keys()

allkeys = sorted(allkeys) # Convert set to sorted list

# Next, try to reduce the list of relations.
# To do this, view the relations obtained in the previous step as elements of the vector space over QQ(g) spanned by graphs found in the relations above.
FQQg = FractionField(QQg)
V = VectorSpace(FQQg,len(allkeys))
vectors = []
for elt in CGs_elements:
    vec = [FQQg(0)]*len(allkeys)
    for k,v in elt.items():
        i = allkeys.index(k)
        vec[i] += FQQg(v)
    vectors.append(V(vec))

W = V.subspace(vectors)
eqns = []
lcm = lambda f,g: (f*g)/g.gcd(f)
for vec in W.basis():
    # simplify the vector by clearing denominators
    denoms = [c.denominator() for c in vec]
    mult = reduce(lcm,denoms,1)
    eqns.append([c*mult for c in vec])

# Output the result.
print("\n\nThe following equations have been found:\n")
for eqn in eqns:
    parts = []
    for i,tup in enumerate(allkeys):
        if eqn[i]!=0:
            parts.append(f"({eqn[i].factor()}) * Gamma{tup}")
    print("  +  ".join(parts) + "\n\n")

print("Here Gamma(r,chi,u) are the r-marked graphs with the following adjacency matrices:\n")
for tup in allkeys:
    print(f"Gamma{tup}:")
    print(seen_graphs[tup[0],tup[1]][tup[2]].matrix())
    print()
