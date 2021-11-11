import numpy as np
from itertools import permutations
from math import factorial
from collections import defaultdict

from .charpoly import charpoly
from .sequences import increasing_injections, integer_sequences, matrices_from_degree_sequence

class MarkedGraph:
    def __init__(self,vertices=None,edges=[],mat=None,marked=0):
        """
        Create a marked graph.
        Keywords:
        marked: number of marked vertices
        mat: adjacency matrix
        
        Instead of mat the following two keywords can be used.
        vertices: number of vertices
        edges: list of tuples representing edges.
        """
        if not mat is None:
            self._mat=np.array(mat,dtype=int)
            self._vertices=len(self._mat)
        elif not vertices is None:
            self._vertices=vertices
            self._mat=np.full((vertices,vertices),0)
            for e in edges:
                self._mat[e[0],e[1]]+=1
            self._mat+=np.transpose(self._mat)
            
        self._marked=marked
        
        # The following variables will be used to cache some invariants of the graph,
        # so they only have to be computed once.
        self._charpoly=None
        self._degree_sequence=None
        self._loop_sequence=None
        
    def degree(self,v):
        """
        Return degree of vertex v
        """
        return sum(self._mat[v,:])
    
    def matrix(self):
        """
        Return adjacency matrix of the graph
        """
        return np.array(self._mat)
    
    def is_contracted(self):
        """
        Returns a boolean that indicates whether the marked graph is contracted.
        """
        return not any(sum(self._mat[v,:])<3 or (sum(self._mat[v,:])==3 and self._mat[v,v]==2) for v in range(self._marked,self._vertices))
    
    def degree_sequence(self):
        """
        Returns a sequence that consists of:
        - first the degrees of the r marked vertices
        - then the degrees of the unmarked vertices, sorted high-to-low
        """
        if self._degree_sequence is None:
            self._degree_sequence = tuple(np.sum(self._mat[:,:self._marked],axis=0)) + tuple(sorted(np.sum(self._mat[:,self._marked:],axis=0), reverse=True))
        return self._degree_sequence
    
    def loop_sequence(self):
        """
        Returns a sequence that consists of:
        - first the numbers of loops at each marked vertex (multiplied by 2)
        - then the numbers of loops at each unmarked vertex, sorted high-to-low (and multiplied by 2)
        """
        if self._loop_sequence is None:
            self._loop_sequence=tuple(self._mat[i,i] for i in range(self._marked)) + tuple(sorted([self._mat[i,i] for i in range(self._marked,self._vertices)],reverse=True))
        return self._loop_sequence
    
    def charpoly(self):
        """
        Returns the coefficients of the characteristic polynomial of the adjacency matrix of the graph
        """
        if self._charpoly is None:
            self._charpoly=charpoly(self._mat)
        return self._charpoly
    
    def marked(self):
        """
        Return number of marked vertices
        """
        return self._marked
    
    def vertices(self):
        """
        Return number of vertices
        """
        return self._vertices
    
    def edges(self):
        """
        Return number of edges
        """
        return np.sum(self._mat)//2
    
    def characteristic(self):
        """
        Return graph characteristic
        """
        return self.vertices()-self.edges()
    
    def is_isomorphic(self,H):
        """
        Returns true if two graphs are isomorphic, false otherwise.
        Algorithm first checks if some basic invariants are equal
        (number of vertices, degree_sequence and loop_sequence,
        marked part and characteristic polynomial of adjacency matrix)
        If yes, then try finding an isomorphism of marked graphs.
        The algorithm used here is not state-of-the-art and more 
        efficient algorithms should exist for comparing graphs.
        For our purposes, however, this is fast enough.
        """
        if self._vertices!=H._vertices or self._marked!=H._marked:
            return False
                
        if self.loop_sequence()!=H.loop_sequence():
            return False
        
        if self.degree_sequence()!=H.degree_sequence():
            return False 
        
        # Check if characteristic polynomials agree:
        if not self.charpoly()==H.charpoly():
            return False
        
        # Check if marked parts are the same:
        if not np.array_equal(self._mat[:self._marked,:self._marked], H._mat[:H._marked,:H._marked]):
            return False
        
        # Otherwise, try to find an isomorphism by permuting rows and columns
        pfix = tuple(range(self._marked))
        for pperm in permutations(range(self._marked,self._vertices)):
            p = pfix+pperm
            if all(self._mat[p[i],p[j]]==H._mat[i,j] for i in range(self._vertices) for j in range(self._vertices)):
                return True
        return False
    
    def __eq__(self,H):
        """
        Returns true if two graphs are isomorphic, false otherwise.
        """
        return self.is_isomorphic(H)
    
    def pushforward(self, phi, s):
        """
        For phi an injective map {1, ..., r} --> {1, ..., s}
        return phi_* of r-marked graph, which is an s-marked graph
        
        s: number of marked vertices of new graph (at least r)
        phi: tuple (phi(1), ..., phi(r)) of length r with coefficients in {1, ..., s}.
        """
        u = self._vertices - self._marked
        newmat = np.zeros((s+u,s+u), dtype=int)
        emap = [_-1 for _ in phi] + list(range(s,s+u))
        newmat[np.ix_(emap,emap)] = self._mat
        return MarkedGraph(marked=s, mat=newmat)
    
    def vertex_is_contracted(self,v):
        """
        Check whether vertex v is contracted.
        """
        if v < self._marked:
            return True
        degv = self.degree(v)
        return not (degv<3 or (degv==3 and self._mat[v,v]==2))
    
    def remove_vertex(self,i):
        """
        Delete vertex i from graph
        """
        if i<self._marked:
            raise Exception('I can\'t remove a marked vertex!')
        else:
            self._mat=np.delete(np.delete(self._mat,i,axis=0),i,axis=1)
            self._vertices-=1
            # reset cache
            self._charpoly=None
            self._loop_sequence=None
            self._degree_sequence=None
            
    def add_edge(self,i,j):
        """
        Add an edge between vertices i and j
        """
        self._mat[i,j]+=1
        self._mat[j,i]+=1
    
    def contract(self,g=None):
        """
        Contract the graph.
        If g is given, returns lambda_g(Gamma)
        """
        if g is not None:
            factor = 1
        while True: # keep looping over all vertices.
            v = next((v for v in range(self._marked,self._vertices) if not self.vertex_is_contracted(v)),None)
            if v is None:
                return factor if (g is not None) else None
            degv = self.degree(v)
            if degv<2:
                self.remove_vertex(v)
                if g is not None and degv==0:
                    factor=0
            else:
                nbs=tuple(w for w in range(self._vertices) if self._mat[v,w]>0)
                if degv==2:
                    if nbs[0]!=v:
                        self.add_edge(nbs[0],nbs[-1])
                    elif g is not None:
                        factor*=(2-2*g)
                    self.remove_vertex(v)
                if degv==3:
                    nb = nbs[0] if nbs[0]!=v else nbs[1]
                    self.add_edge(nb,nb)
                    self.remove_vertex(v)


_CGposcache = dict()
def CGpos(r,chi,u):
    """
    Returns a full list of (graphs representing) the equivalence classes
    of r-marked graphs of characteristic chi and u unmarked vertices,
    for which each marked vertex has positive degree.
    
    This function is used as a helper function for CG(r,chi,u).
    """
    # Check if we computed this before
    if (r,chi,u) in _CGposcache:
        return _CGposcache[r,chi,u]
    
    L = list()
    # First, we generate all marked graphs whose marked vertices
    # all have positive degree.
    e = r + u - chi # number of edges
    # Loop over S1: sum of degrees of marked vertices
    for S1 in range(2*e-3*u+1): 
        S2 = 2*e - S1 # sum of degrees of unmarked vertices
        for deg_seq_1 in integer_sequences(L=r,S=S1,m=(1,)*r):
            for deg_seq_2 in integer_sequences(L=u,S=S2,nondecr=True,m=(3,)*u):
                deg_seq = deg_seq_1 + deg_seq_2
                for M in matrices_from_degree_sequence(deg_seq):
                    G = MarkedGraph(marked=r, mat=M)
                    if G.is_contracted() and not G in L:
                        L.append(G)
    _CGposcache[r,chi,u] = L
    return L

_CGcache = dict()                    
def CG(r,chi,u=None):
    """
    Returns a full list of (graphs representing) the equivalence classes
    of r-marked graphs of characteristic chi and u unmarked vertices.
    """
    # If u is not set, return CG(r,chi)
    if u is None:
        ret = []
        for u in range(2*(r-chi)+1):
            ret += CG(r,chi,u)
        return ret
    
    # Check if we computed this before
    if (r,chi,u) in _CGcache:
        return _CGcache[r,chi,u]
    
    L = list()
    # Loop over the number of marked vertices with positive degree.
    for rpos in range(0,r+1):
        Lnew = list()
        # Loop over all injective maps {1, ..., rpos} --> {1, ..., r}
        for phi in increasing_injections(rpos,r):
            for G0 in CGpos(rpos,chi+rpos-r,u):
                G = G0.pushforward(phi,r)
                if G not in Lnew:
                    Lnew.append(G)  
        L += Lnew
    
    _CGcache[r,chi,u] = L
    return L

