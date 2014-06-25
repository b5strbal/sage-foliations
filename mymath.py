"""
General supporting math functions.

AUTHORS:

- Balazs Strenner (2014-06-16): initial version

EXAMPLES::


"""

#*****************************************************************************
#       Copyright (C) 2014 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from collections import deque
from sage.matrix.constructor import matrix, vector
from sage.symbolic.ring import var
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import \
    PolynomialRing
from constants import epsilon
from sage.functions.other import floor 

def mod_one(x):
    """
    Return ``x`` mod 1.

    INPUT:

    - ``x`` -- a real number

    OUTPUT:

    ``x`` mod 1, a number in `[0,1)`.

    TESTS::

        sage: from sage.dynamics.foliations.mymath import mod_one
        sage: mod_one(2.5)
        0.500000000000000
        sage: mod_one(-1.7)
        0.300000000000000
        sage: mod_one(7/6)
        1/6
        sage: mod_one(-1/6)
        5/6
        sage: a = QQbar(sqrt(2)); a
        1.414213562373095?
        sage: mod_one(a)
        0.4142135623730951?

    """
    return x - floor(x)

class OrientedGraph(SageObject):
    r"""
    Simple digraph class for checking primitivity of a matrix.

    Vertices of the graph are represented by numbers between 0 and
    `n-1` where `n` is the number of vertices.

    INPUT:
        
    - ``square_matrix`` -- a square matrix with non-negative integer
    entries, the adjacency matrix of the graph.

    """
    def __init__(self, square_matrix):
        r"""
        Constructor.
        
        """
        m = square_matrix
        assert(m.is_square())
        for i in range(m.nrows()):
            for j in range(m.ncols()):
                asssert(m[i,j] >= 0)
                
        self._n = m.nrows()
        self._adj_list = matrix_to_adj_list(m)
        self._backward_adj_list = matrix_to_adj_list(m.transpose())

    def bfs_distances(self, v, backward = False):
        """
        Return the distances of vertices measured from ``v``.

        The distance from itself, 0, is also included.
        
        INPUT:
        
        - ``v`` -- an integer, representing a vertex

        - ``backward`` -- (default: False) if False, the distances are
          measured FROM ``v``, otherwise TO ``v``.

        OUTPUT:

        the list of the distances (non-negative integers)

        """
        adj_list = self._backward_adj_list if backward else self._adj_list
        done = set([v])
        queue = deque([v])
        distances = [None] * self._n
        distances[v] = 0
        while len(queue) > 0:
            u = queue.popleft()
            for w in adj_list[u]:
                if w not in done:
                    queue.append(w)
                    done.add(w)
                    distances[w] = distances[u] + 1
        return distances

    def is_primitive(self):
        """Decide if the adjacency matrix is primitive.

        OUTPUT:

        True of False, depending on whether the matrix is primitive or
        not

        ALGORITHM:

        `Leegard
        <http://ir.library.oregonstate.edu/xmlui/bitstream/handle/1957/31568/LeegardAmandaD2003.pdf?sequence=1>`_
        [LEE2002]_ has given a fast algorithm for determining
        primitivity. The idea is that it is sufficient to calculate
        the lengths of certain cycles containing a fixed vertex ``v``,
        and the matrix adjacency matrix is primitive if and only if
        the gcd of these lengths is 1.

        REFERENCES:

        .. [LEE2002] Amanda D. Leegard. A Fast Algorithm for
        Determining Primitivity of an `n\times n` non-negative
        matrix. Thesis, 2002.

        """
        P = self.bfs_distances(0, backward = True)
        Q = self.bfs_distances(0, backward = False)
        if None in P or None in Q: # not strongly connected
            return False
        C = {Q[i] + 1 + P[j] for i in xrange(self._n)
             for j in self._adj_list[i]} # cycle lengths including 0
        # checking for self-loops
        if any(i in self._adj_list[i] for i in xrange(self._n)):
            return True
        # once it is strongly connected, it is primitive if and only if
        # the gcd of cycle lengths is 1
        return gcd(list(C)) == 1


def fixed_charpoly(M, variables):
    """
    Calcualate the characterictic polynomials of a matrix with entries
    in a group ring.

    The problem is that Sage doesn't do a good job at calculating
    characteristic polynomials of matrices over Symbolic Ring in
    certain cases. To overcome this, we create a new matrix over a
    Polynomial Ring where inverses of the variables are replaced with
    new variable, calculate the charpoly, then substitute back.

    INPUT:

    - ``M`` -- a square matrix over Symbolic Ring

    - ``variables`` -- the list of variables appearing in the entries
      of ``M``

    OUTPUT:

    the characterictic polynomial of ``M``

    """
    sub_vars = {v : var(str(v) + 'inv') for v in variables}
    # print sub_vars
    ring = PolynomialRing(ZZ, variables + sub_vars.values()) \
           if len(variables) > 0 else ZZ
    new_M = matrix(ring, M.nrows(), M.ncols())
    for i in xrange(M.nrows()):
        for j in xrange(M.ncols()):
            exp = M[i,j]
            a = exp.numerator()
            b = exp.denominator()
            # print a, b, b.subs(sub_vars)
            new_M[i,j] = a * b.subs(sub_vars)
    # print new_M
    # print new_M.parent()
    poly = new_M.charpoly('x')
    back_sub = {str(sub_vars[v]) : 1/v for v in sub_vars}
    # print poly
    # print poly.parent()
    # print back_sub
    poly = poly.subs(**back_sub)
    return poly

            
def matrix_to_adj_list(square_matrix):
    n = square_matrix.nrows()
    return [[j for j in xrange(n) if square_matrix[i,j] > 0]
            for i in xrange(n)]
        


def is_perron_frobenius(square_matrix):
    # print square_matrix
    # print square_matrix**50
    g = OrientedGraph(square_matrix)
    # print g.is_primitive()
    return g.is_primitive()

class NoPFEigenvectorError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

def pf_eigen_data(square_matrix, field):
    m = matrix(square_matrix, field)
    evr = m.eigenvectors_right()
    largest = max(evr, key = lambda x: abs(x[0]))
    eigenvalue = largest[0]
    if abs(eigenvalue.imag()) > epsilon or eigenvalue.real() < 1 + epsilon:
        raise NoPFEigenvectorError("Dominant eigenvalue is not a real number"
                                   " greater than 1.")
    if largest[2] > 1:
        raise NoPFEigenvectorError("Dominant eigenvalue has multiplicity "
                                   "larger than one.")
    v = make_positive(largest[1][0])  # eigenvector
    if v == -1: # cannot made positive
        raise NoPFEigenvectorError("Dominant eigenvector doesn't lie in the "
                                   "positive quadrant.")
    v /= sum(v)
    return (eigenvalue, v)
    



def make_positive(vec):
    """
    Returns a parallel vector to a given vector with all positive
    coordinates if possible.

    INPUT:

    - ``vec`` - a vector or list or tuple of numbers

    OUTPUT:

    - vector - a vector with positive coordinates if this is possible,
      otherwise -1

    EXAMPLES:

    If all coordinates of the input vector are already positive, the
    same vector is returned::

        sage: from sage.dynamics.foliations.foliation import make_positive
        sage: v = vector([1, 3, 5])
        sage: make_positive(v) == v
        True

    If all are negative, its negative is returned::

        sage: make_positive((-1, -5))
        (1, 5)

    Even if the coordinates are complex, but have very small imaginary
    part as a result of an approximate eigenvector calculation for 
    example, the coordinates are treated as real::

        sage: make_positive((40.24 - 5.64e-16*I, 1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)
        sage: make_positive((-40.24 - 5.64e-16*I, -1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)

    If there is a complex coordinate which is not negligible, -1 is
    returned::

        sage: make_positive((-40.24 - 5.64e-6*I, -1.2 + 4.3e-14*I))
        -1

    If one coordinate is zero, or very close to zero, -1 is returned::

        sage: make_positive((1, 0, 2))
        -1
        sage: make_positive((-40.24e-15*I - 5.64e-16*I, -1.2))
        -1

    If there are both negative and positive coordinates, -1 is 
    returned::

        sage: make_positive((-3, 4, 5))
        -1
        
    """
    if any(abs(x) < epsilon for x in vec) or \
        any(abs(x.imag()) > epsilon for x in vec):
            return -1
    newvec = vector([x.real() for x in vec])
    if vec[0].real() < 0:
        newvec = -newvec
    if any(x < 0 for x in newvec):
        return -1
    return newvec

# def get_good_eigendata(transition_matrix, is_twisted):
#     """
    
#         sage: from sage.dynamics.foliations.foliation import get_good_eigendata


#     """
#     ev = transition_matrix.eigenvectors_right()
#     ret_list = []
#     for x in ev:
#         if abs(x[0].imag()) < epsilon and x[0].real() > 0 \
#                 and abs(x[0].real() - 1) > epsilon:
#             for vec in x[1]:
#                 newvec = make_positive(vec)
#                 if newvec != -1:
#                     norm = sum([abs(y) for y in newvec])
#                     if is_twisted:
#                         norm -= abs(newvec[-1])
#                     normalized_vec = [abs(y / norm) for y in newvec]
#                     ret_list.append((x[0].real(), normalized_vec))
#     return ret_list





