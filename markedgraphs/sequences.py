import numpy as np

def integer_sequences(L, S, nondecr=False, m=None, M=None):
    """
    Generate sequences of non-negative integers.
    
    Parameters:
    L:         the length of the sequences
    S:         the sum of the integers in each sequence
    
    Optional parameters:
    nondecr:   (boolean) return only non-decreasing sequences (default: False)
    m:         tuple of length L; gives lower bounds for coefficients of list (default: None)
    M:         tuple of length L; gives upper bounds for coefficients of list (default: None)
    """
    
    # If M and m are not given, use the following defaults.
    if M is None:
        M = (S,)*L
    if m is None:
        m = (0,)*L
        
    # If length=0 and sum=0 then yield the empty tuple.
    # Otherwise, yield nothing.
    if L==0:
        if S==0:
            yield tuple()
    # If length=1 and sum lies within given boundaries, yield 1-tuple.
    elif L==1:
        if m[0]<=S<=M[0]: 
            yield (S,)
    # If length>1, then loop through possible values for first coefficient,
    # and recursively generate (L-1)-tuples that give the remaining coefficients.
    elif L>1:
        for first in range(m[0], min(S,M[0])+1):
            m_next = m[1:] if nondecr==False else (first,)+m[2:]
            for tail in integer_sequences(L=L-1, S=S-first, nondecr=nondecr, m=m_next, M=M[1:]):
                yield (first,)+tail

def matrices_from_degree_sequence(deg_seq=tuple()):
    """
    Generate all matrices that can occur as adjacency matrices
    of r-marked graphs whose degree_sequence() equals deg_seq.
    
    These are symmetric (r×r)-matrices with non-negative integer coefficients
    with even coefficients on the diagonal, whose row sum equals deg_seq.
    """
    L = len(deg_seq)
    
    M = np.zeros((L,L), dtype=int)
    if L==0:
        yield M
    elif L>0:
        # Range over possible values for top left entries:
        for top_left in range(0,deg_seq[0]+1,2):
            M[0,0] = top_left
            # Range over sequences that make up the rest of the first row:
            for top_row_rest in integer_sequences(L=L-1, S=deg_seq[0]-top_left, M=deg_seq[1:]):
                M[0,1:] = M[1:,0] = top_row_rest
                # Compute the row sum of the remaining bottom right square
                row_sum_rem = tuple(deg_seq[i] - M[0,i] for i in range(1,L))
                
                # Loop over all possible bottom right squares:
                for BRS in matrices_from_degree_sequence(deg_seq=row_sum_rem):
                    M[1:,1:]=BRS
                    yield M.copy()

def increasing_injections(r,s,m=1):
    """
    Generates all increasing tuples of length r
    with values in {1, ..., s}.
    
    Optional argument:
    m: minimum value of first entry (default=1, used for recursion)
    """
    if r==0:
        yield tuple()
    elif r==1:
        for i in range(m,s+1):
            yield (i,)
    else:
        for first in range(m,s-r+2):
            for tail in increasing_injections(r-1,s,m=first+1):
                yield (first,)+tail
                    