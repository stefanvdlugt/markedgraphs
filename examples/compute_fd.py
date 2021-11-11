# The following script can be used to compute the polynomials f_{d,u} and f_d, where:
# f_{d,u}(r) is the number of contracted r-marked graphs of characteristic r-d with u unmarked vertices, and
# f_d(r) is the number of contracted r-marked graphs of characteristic r-d.
# Requires sympy.

# allow import of markedgraphs module from parent directory:
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from sympy import latex, symbols, prod, expand
from markedgraphs import CG

d = int(input("Enter d: "))
ui = input("Enter u (leave empty to compute f_{d,u} for all 0<=u<=2d): ")

# ulist contains all values of u for which we should compute f_{d,u}.
ulist = [int(ui)] if ui else range(2*d+1)

# fdu[u] equals f_{d,u}.
# firstvalues[u] lists the values f_{d,u}(r) for 0<=r<=2d-u.
fdu = dict()
firstvalues = dict()

rr = symbols('r')

for u in ulist:
    firstvalues[u] = [len(CG(r,r-d,u)) for r in range(2*d-u+1)]
    
    # Compute f_{d,u} using Lagrange interpolation:
    lagrange_basis = []
    for i in range(2*d-u+1):
        li = prod((rr-j)/(i-j) for j in range(2*d-u+1) if j!=i)
        lagrange_basis.append(li)
    fdu[u] = sum(firstvalues[u][i]*lagrange_basis[i] for i in range(2*d-u+1))
    
    

# Output results:
print()
for u in ulist:
    print(f"First values of f_{{{d},{u}}}:",firstvalues[u])


print("\n\nLagrange interpolation yields:")
for u in ulist:
    print(f"f_{{{d},{u}}} =", latex(expand(fdu[u])))
    
# If u was not given, also output f_d
if not ui:
    print()
    fd = sum(fdu[u] for u in ulist)
    print(f"f_{{{d}}} =", latex(expand(fd)))
