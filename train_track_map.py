r"""


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


from base import *
from collections import namedtuple
from sage.matrix.constructor import matrix, vector
from sage.rings.integer_ring import ZZ
from train_track import TrainTrackPath

class TrainTrackMap(namedtuple("TrainTrackMap", "domain, codomain,"
                     "vertex_map, edge_map, coding_list")):

    def __repr__(self):
        return "Train track map from " + repr(self.domain) + " to " +\
            repr(self.codomain)

    def __eq__(self, other):
        return isinstance(other, TrainTrackMap) and \
            self.domain == other.domain and \
            self.codomain == other.codomain and \
            self.vertex_map == other.vertex_map and \
            self.edge_map == other.edge_map

    def __mul__(self, other):
        # f*g means g comes first, then f
        # print self.domain.foliation
        # print self.codomain.foliation
        # print other.domain.foliation
        # print other.codomain.foliation
        # print self.edge_map
        # print other.edge_map
        new_vertex_map = {v:self.vertex_map[other.vertex_map[v]]
                          for v in other.vertex_map}
        new_edge_map = {e:self._map_path(other.edge_map[e])
                        for e in other.edge_map}

        new_map = TrainTrackMap(domain = other.domain,
                                codomain = self.codomain,
                                vertex_map = new_vertex_map,
                                edge_map = new_edge_map,
                                coding_list = other.coding_list + self.coding_list)
        # if new_map.edge_matrix() !=
        #        self.edge_matrix() * other.edge_matrix():
        # print new_map.edge_matrix(), '\n'
        # print other.edge_matrix() * self.edge_matrix(), '\n'
        # print self.edge_matrix(), '\n'
        # print other.edge_matrix(), '\n'
        # print other.domain._index_to_edge
        # print other.codomain._index_to_edge
        assert(new_map.edge_matrix() ==
               other.edge_matrix() * self.edge_matrix())
        # assert(new_map.small_matrix() ==
        #        other.small_matrix() * self.small_matrix())

        return new_map

    def _map_path(self, path):
        new_path = []
        for oe in path:
            p = self.edge_map[oe.edge]
            if oe.direction == 1:
                new_path.extend(p)
            else:
                new_path.extend(p.reversed())
        return TrainTrackPath(new_path)

    def _cohomology_kernel(self):
        # this only works if domain and codomain are combinatorically same
        # the underlying foliations need not be the same
        tt = self.domain

        m = self.edge_matrix(SIGNED)
        # m = matrix([tt.path_to_vector(self.edge_map[edge], SIGNED)
        #             for edge in tt.edges()])
        m = m.transpose() - matrix.identity(ZZ, m.nrows())
        mlist = m.rows()
        mlist.extend(tt.kernel_from_singularities())
        return matrix(mlist).right_kernel_matrix()


    def is_self_map(self):
        return self.domain == self.codomain

    def action_on_cohomology(self, is_punctured):
        # when the surface is nonorientable, this matrix should be
        # shrinking with eigenvalue 1/alpha, since it is a pullback.
        basis = self.codomain.basis_of_cohomology(is_punctured)
        new_vectors = basis * self.edge_matrix(SIGNED).transpose()

        new_vectors *= matrix(ZZ, self.domain.cochain_matrix().inverse())
        # the rows represent the range of the linear map, so we transpose it
        return new_vectors.matrix_from_columns(range(new_vectors.nrows())).transpose()

    def invariant_cohomology(self, is_punctured):
        assert(self.is_self_map())
        action = self.action_on_cohomology(is_punctured)
        m = action - matrix.identity(ZZ, action.nrows())
        kernel = m.right_kernel_matrix()
        basis = self.domain.basis_of_cohomology(is_punctured)

        result = kernel * basis
        return result



    def small_matrix(self):
        # The domain and the codomain should be both one-sided or both
        # two-sided, otherwise the matrix won't be a square matrix
        # Usually this is not a problem, since we only call this method
        # is the underlying permutations are the same which is a much stronger
        # condition.
        m = self.domain.matrix_to_reduce_dimension()
        # print m
        # print self.domain.foliation
        # print self._edge_matrix, '\n'
        result = m.transpose() * self.edge_matrix()
        # print result, '\n'
        result = result.matrix_from_columns(range(
            self.codomain.matrix_to_reduce_dimension().ncols()))
        # print result, '\n'
        
        return result

    def edge_matrix(self, is_signed = UNSIGNED):
        if not hasattr(self, '_edge_matrix'):
            self._edge_matrix = [None, None]

        if self._edge_matrix[is_signed] == None:
            self._edge_matrix[is_signed] = matrix([self.codomain.path_to_vector(
                self.edge_map[edge], is_signed)
                                                   for edge in self.domain.edges()])
            assert(set(self.edge_map.keys()) == set(self.domain.edges())) 
        # if abs(self.small_matrix().det()) != 1:
        #     print self._edge_matrix
        #     print self.small_matrix().det()
        #     print self
        #     print self.small_matrix()
        #     print self.small_matrix().charpoly().factor()
        # assert(abs(self.small_matrix().det()) == 1)
        return self._edge_matrix[is_signed]


    def get_PF_weights(self, field):
        return pf_eigen_data(self.edge_matrix().transpose(), field)
    
        # if not self.is_self_map():
        #     raise ValueError("The train track map should be a self-map.")
        # m = tt_map.edge_matrix().transpose()
        # try:
        #     return pf_eigen_data(m, RDF)
        # except NoPFEigenvectorError as ex:
        #     # print "No PF Eigendata: ", ex
        #     raise ValueError("The map is not pseudo-Anosov.")

    def get_PF_interval_lengths(self, field):
        return list(pf_eigen_data(self.small_matrix().transpose(), field)[1])        

        # m = self.small_matrix().transpose()
        # try:
        #     return list(pf_eigen_data(m, RDF)[1])
        # except NoPFEigenvectorError as ex:
        #     # print "No PF Eigendata: ", ex
        #     continue

        
    def is_pseudo_anosov(self):
        return is_perron_frobenius(self.edge_matrix())






def tt_map_from_codings(start_tt,coding_list):
    from train_track_map import TrainTrackMap
    from foliation import SymmetryCoding
    from transverse_curve import RestrictionCoding, TransverseCurve, Coding
    tt_map_so_far = start_tt.identity_map()
    for coding in coding_list:
        if isinstance(coding, SymmetryCoding):
            tt_map = tt_map_so_far.domain.transform(*coding)
        elif isinstance(coding, RestrictionCoding):
            tc = TransverseCurve(coding.foliation, coding.coding)
            fol, tt_map = tc.new_foliation()
        elif isinstance(coding, Coding):
            tc = TransverseCurve(fol, coding)
            fol, tt_map = tc.new_foliation()
        else:
            assert(False)

        tt_map_so_far = tt_map_so_far * tt_map            

    return tt_map_so_far




from base import MyException
class NoPFEigenvectorError(MyException):
    pass

def pf_eigen_data(square_matrix, field):
    """Return the Perron-Frobenius eigenvalue and -vector of a matrix.

    A PF eigenvalue is an eigenvalue which is a positive real number,
    greater in absolute value than all the other eigenvalues, and
    greater than 1 as well. It is of multiplicity one, and the
    coordinates of the corresponding eigenvector are all positive.
    
    INPUT:

    - ``square_matrix`` -- a square matrix with integer entries (not
      necessarily positive).

    - ``field`` -- the field in which the eigendata is
      calculated. Whether it is ``QQbar`` of ``RDF``, the return value
      will contain exact algebraic numbers in ``QQbar``, or standard
      floating point numbers. The latter choice is much faster.

    OUTPUT:

    a tuple ``(a,v)``, where ``a`` is the PF eigenvalue and ``v`` is
    the corresponding eigenvector. If no such eigenvalue exists,
    ``NoPFEigenvectorError`` is raised.

    """

    m = matrix(square_matrix, field)
    evr = m.eigenvectors_right()
    largest = max(evr, key = lambda x: abs(x[0]))
    eigenvalue = largest[0]
    if abs(eigenvalue.imag()) > EPSILON or eigenvalue.real() < 1 + EPSILON:
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
    if any(abs(x) < EPSILON for x in vec) or \
        any(abs(x.imag()) > EPSILON for x in vec):
            return -1
    newvec = vector([x.real() for x in vec])
    if vec[0].real() < 0:
        newvec = -newvec
    if any(x < 0 for x in newvec):
        return -1
    return newvec








from collections import deque
# from sage.structure.sage_object import SageObject

class OrientedGraph():
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

        A matrix being primitive is the same as being Perron-Frobenius.

        OUTPUT:

        True of False, depending on whether the matrix is primitive or
        not.

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

            
def matrix_to_adj_list(adj_matrix):
    """Create an adjacency list from an adjacency matrix.

    INPUT:

    - ``adj_matrix`` -- a square matrix with non-negative integer
      entries, imagined and the adjacency matrix of a graph.

    OUTPUT:

    The adjacenct list as a list of lists. The i'th list in the list
    contains the vertices to which there is an oriented edge from the
    i'th vertex. Multiple edges from the same pair vertices are
    ignored, but this is not okay, because this is used only for
    deciding the primitivity of the adjacency matrix.
    
    """

    n = adj_matrix.nrows()
    return [[j for j in xrange(n) if adj_matrix[i,j] > 0]
            for i in xrange(n)]
        


def is_perron_frobenius(square_matrix):
    """Decides if a square matrix is Perron-Frobenius.

    A non-negative integer square matrix is called Perron-Frobenius if
    some power of it contains only positive entries.

    INPUT:

    - ``square_matrix`` -- a square matrix with non-negative integer
      entries

    OUTPUT:

    true of false accoring to whether the matrix is Perron-Frobenius

    """
    g = OrientedGraph(square_matrix)
    return g.is_primitive()



