# Input: integers r, s, G
# Constructs elements in the kernel of the map R[g][CG(r)] --> R^*(C_g^r)
# and sends these elements through the map R[g][CG(r)] --> R[g][G(s)] --> R[g][CG(s)]
# to obtain elements in the kernel of the map R[g][CG(s)] --> R^*(C_g^s).

from sage.all import CombinatorialFreeModule, PolynomialRing, QQ

# allow import of markedgraphs module from parent directory:
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from markedgraphs import CG, MarkedGraph
from collections import defaultdict
from math import gcd
from functools import reduce

if len(sys.argv)==4:
    r = int(sys.argv[1])
    s = int(sys.argv[2])
    G = int(sys.argv[3])
else:
    r = int(input('Enter r: '))
    s = int(input('Enter s: '))
    G = int(input('Enter G: '))

d = 2*(G+1-r+s) # degree of resulting forms on C_g^s

def lcm(a,b):
    return a*b//gcd(a,b)

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
    
def prim_factor(L):
    """
    Take a vector of rationals and return a rational q
    such that multiplying the vector with q yields
    a primitive vector of integers.
    """
    if len(L)==0:
        return []
    denoms = [(_).denominator() for _ in L]
    mult = reduce(lcm,denoms)
    L = [int(_*mult) for _ in L]
    div = reduce(gcd,L)
    return mult/div
    
    
# Construct the polynomial ring Q[g]
QQg = PolynomialRing(QQ,'g')
g = QQg.gens_dict()['g']

# Construct the polynomial ring Z[g][Gamma_{i,j}] with 1 <= i <= j <= r.
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
relations = (Wr**(G+1)).coefficients()

# For each of the coefficients found, take the image 
# under the homomorphism Z[g][Gamma_{i,j}] --> Z[g][CG(r)] --> Z[g][G(s)] --> Z[g][CG(s)]
# We'll represent elements of Z[g][CG(s)] as dicts. 
# The keys correspond to elements of CG(s), and its values are the coefficients in Z[g].

# edge[k] = Gamma_{i,j} iff Gamma_{i,j} is the k-th generator of R.
edge = {R.gens().index(Gamma[i,j]): (i,j) for i in range(r) for j in range(i,r)}
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
    #elt = {k:v for k,v in elt.items() if v!=0}
    if elt:
        CGs_elements.append(elt)

        
# We have now found some vectors;
# we'll let sage put them in echelon forms
# so we have no linearly dependent ones.
QQg = PolynomialRing(QQ,'g')
all_tups = set()
for elt in CGs_elements:
    all_tups |= set(elt.keys())
genlist = [f'Gamma_{a}_{b}_{c}' for a,b,c in all_tups]
MM = CombinatorialFreeModule(QQg, genlist)
vectors = list()
for elt in CGs_elements:
    vector = MM(0)
    for (a,b,c),v in elt.items():
        vector += MM(f'Gamma_{a}_{b}_{c}')*QQg(v)
    vectors.append(vector)

# Output results
print( "The following elements lie in the kernel of the homomorphism\n")
print(f"   alpha_{{{s}}}: Z[CG({s})] --> R^{{{d}}}(C_g^{{{s}}})\n")
print(f"for all 2 <= g <= {G}:\n")
    
printed_graphs = set()
for vector in MM.echelon_form(vectors):
    parts = []
    mult = prim_factor([cc for c in vector.coefficients() for cc in c.coefficients()])
    vector = mult*vector
    for a,b,c in sorted(all_tups):
        vc = vector.coefficient(f'Gamma_{a}_{b}_{c}')
        if vc!=0:
            printed_graphs.add((a,b,c))
            parts.append(f'( {QQg(vc).factor()} * Gamma({a},{b},{c}) )')
    print('  +  '.join(parts), end='\n\n')

print("Here Gamma(r,chi,u) are the r-marked graphs with the following adjacency matrices:\n")
for tup in sorted(printed_graphs):
    print(f"Gamma{tup}:")
    print(seen_graphs[tup[0],tup[1]][tup[2]].matrix())
    print()

