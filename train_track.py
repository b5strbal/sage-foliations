from sage.structure.sage_object import SageObject
from mymath import mod_one
from collections import namedtuple
from sage.matrix.constructor import matrix, vector
from sage.rings.integer_ring import ZZ
from constants import *
from sage.misc.misc_c import prod


class OrientedEdge(namedtuple('OrientedEdge', 'edge, direction')):
    def __new__(cls, edge, direction):
        return super(OrientedEdge, cls).__new__(cls, edge, direction)
        
    def reversed(self):
        return OrientedEdge(self.edge, -self.direction)

    # def endpoint(self, end):
    #     if (end == 1) == (self.direction == 1):
    #         return self.edge[0]
    #     return self.edge[1]

    def start(self):
        if self.direction == 1:
            return self.edge[0]
        return self.edge[1]

    def end(self):
        if self.direction == 1:
            return self.edge[1]
        return self.edge[0]


class TrainTrack(SageObject):
    def __init__(self, foliation):
        self._sample_fol = foliation

        self._from = {'pair' : {}, 'center' : {}}
        self._vertex_to_index = {}
        self._index_to_vertex = foliation.intervals()
        self._edge_to_index = {}
        self._index_to_edge = []
        vertex_count = 0
        edge_count = 0
        for interval in foliation.intervals():
            self._vertex_to_index[interval] = vertex_count
            vertex_count += 1
            if not interval in self._from['pair'].values():
                pair = interval.pair(foliation)
                self._from['pair'][interval] = pair
                edge = (interval, pair, 'pair')
                self._edge_to_index[edge] = edge_count
                self._index_to_edge.append(edge)
                edge_count += 1

        for interval in foliation.intervals():
            side = interval.endpoint_side(LEFT, foliation)
            containing_int = foliation.in_which_interval(\
                                    interval.endpoint(LEFT, foliation), (side + 1) % 2)
            self._from['center'][interval] = containing_int
            edge = (interval, containing_int, 'center')
            self._edge_to_index[edge] = edge_count
            self._index_to_edge.append(edge)
            edge_count += 1

        self._coboundary_map = [None, None]
        self._tree_edge_indices = None
        
        self._degrees = [0] * len(self._index_to_vertex)
        for edge in self.edges():
            for i in [0,1]:
                self._degrees[self._vertex_to_index[edge[i]]] += 1
                        
        self._degree_product = prod(self._degrees)

        n = foliation.num_intervals(TOP)
        d1, d2 = self._degrees[:n], self._degrees[n:]
        self._checksum = sum([abs(d[i]-d[(i+1)%len(d)]) for d in [d1, d2]
             for i in xrange(len(d))])

        self._invariants = (self._degree_product, self._checksum, n)

    def __eq__(self, other):
        # print self._index_to_edge
        # print other._index_to_edge
        # print self._index_to_edge == other._index_to_edge
        # print '\n'
        
        return isinstance(other, TrainTrack) and \
            self.invariants() == other.invariants() and \
            self._index_to_edge == other._index_to_edge and \
            self.sample_fol().permutation() == other.sample_fol().permutation()

        # The flips of the sample foliation should be the same as well,
        # otherwise the underlying surfaces might not be homeomorphic.


    def _repr_(self):
        s = "Train track with degrees "
        s += repr(self._degrees) + " and invariants "
        s += repr(self.invariants())
        return s


    def identity_map(self):
        from train_track_map import TrainTrackMap
        return TrainTrackMap(self, self,
                             {v:v for v in self.vertices()},
                             {e:TrainTrackPath.from_edge(e) for e in self.edges()},
                             coding_list = [])


    def transform(self, interval, hdir):
        return self.sample_fol().transform(interval, hdir)[1]

    def invariants(self):
        return self._invariants

    def sample_fol(self):
        return self._sample_fol

    # def degree_product(self):
    #     return self._degree_product

    # def checksum(self):
    #     return self._checksum

    def vertices(self):
        return self._index_to_vertex

    def vertex_to_index(self, vertex):
        return self._vertex_to_index[vertex]


    def edges(self):
        return self._index_to_edge # list

    def kernel_from_singularities(self):
        circles = self.sample_fol().paths_around_singularities
        return matrix([self.path_to_vector(circle, SIGNED) for circle in circles])

        # the row vectors might not be primitive, so to scale them down
        # we echelonize
        # return m.echelon_form(include_zero_rows = False)

    def coboundary_map_from_vertices(self, cochain_type):
        if self._coboundary_map[cochain_type] == None:             
            m = matrix(ZZ, len(self._index_to_edge), len(self._index_to_vertex))
            if cochain_type == COHOMOLOGY:
                for i in range(len(self._index_to_edge)):
                    e = self._index_to_edge[i]
                    m[i, self._vertex_to_index[e[0]]] -= 1
                    m[i, self._vertex_to_index[e[1]]] += 1
            elif cochain_type == TRAIN_TRACK_MODULE:
                for i in range(len(self._index_to_edge)):
                    e = self._index_to_edge[i]
                    if e[2] == 'center':
                        x = 1
                    else:
                        x = -1
                    m[i, self._vertex_to_index[e[0]]] += x
                    m[i, self._vertex_to_index[e[1]]] += x
            self._coboundary_map[cochain_type] = m

        return self._coboundary_map[cochain_type]

    def _compute_edge_cocycles(self):
        C = self.coboundary_map_from_vertices(COHOMOLOGY)
        D, U, V = C.smith_form() # D = U * C * V

        # Modifying U in the smith form such that D is a diagonal matrix
        # where the non-zero rows are in the end instead of at the beginning
        U = matrix(list(reversed(U.rows())))

        # make a matrix such that the rows generate the edge cocycles
        # the last num_vertices-1 rows are coboundaries
        U = matrix(ZZ, U.transpose().inverse())

        # separating the matrix to coboundaries and its complement 
        # and echelonize the coboundaries
        k = len(self.vertices()) - 1
        coboundaries = U[-k:].echelon_form()
        basis = U[:-k]

        # adjusting the basis such that all entries are zero in the columns
        # definied by the pivots of coboundaries (which in turn define the
        # edges of the minimal subtree if the train track where all the
        # cohomology representatives vanish.
        pivots = coboundaries.pivots()
        X = basis.matrix_from_columns(pivots)
        basis -= X * coboundaries

        # now we the subspace generated by the rows of basis into two
        # the subspace vanishing under the coboundary map C^1 -> C^2
        # and its complement
        m = basis * self.kernel_from_singularities().transpose()

        D, U, V = m.smith_form()

        basis = U * basis
        basis = matrix(list(reversed(basis.rows())))
        
        self._cochains = matrix(basis.rows() + coboundaries.rows())
        self._basis_of_cohomology_punctured = basis
        self._basis_of_cohomology_closed = basis[:-m.rank()]
        
        
    def basis_of_cohomology(self, is_punctured):
        if not hasattr(self, "_cochains"):
            self._compute_edge_cocycles()
        if is_punctured:
            return self._basis_of_cohomology_punctured
        return self._basis_of_cohomology_closed

    def cochain_matrix(self):
        if not hasattr(self, "_cochains"):
            self._compute_edge_cocycles()
        return self._cochains


                    
    def matrix_to_reduce_dimension(self):

        # BUG: when the twist is longer than the first interval, this, and the 
        # small matrix calculation is buggy
        
        if hasattr(self, '_reducing_matrix'):
            return self._reducing_matrix
        from copy import copy
        from sage.rings.rational import Rational
        X = self.coboundary_map_from_vertices(TRAIN_TRACK_MODULE)
        # print 'start'
        # print X, '\n'
        # flipping the matrix over vertically (so the echelonizing works
        # in the way we want) and transpose
        X = matrix(list(reversed(X.rows()))).transpose()
        X = X.echelon_form()
        # print X, '\n'
        
        # if the train track is non-orientable, the last row will be a
        # multiple of two, so we have to normalize it and re-echelonize.
        n = X.nrows()
        X = copy(X).with_row_set_to_multiple_of_row(n-1, n-1, Rational('1/2'))
        X = X.echelon_form()
        # print X, '\n'

        # now we need to cut off the identity part of the matrix and
        # append a small block of identity in the different direction
        identity_size = X.nrows() if self.sample_fol().is_bottom_side_moebius() \
                        else X.nrows() - 1
        assert((X[X.nrows()-1, X.nrows()-1] == 1) == (identity_size == X.nrows()))
        X = matrix(list(reversed(X.rows()[:identity_size])))
        X = -matrix(list(reversed(X.columns()[identity_size:]))).transpose()
        X = matrix(matrix.identity(X.ncols()).rows() + X.rows())
        # print X, '\n'
        self._reducing_matrix = X
        return X

    def get_edge_from(self, from_vertex, edge_type):
        return (from_vertex, self._from[edge_type][from_vertex], edge_type)
        
    def get_center_edge(self, from_vertex, end):
        if end == LEFT:
            return OrientedEdge(self.get_edge_from(from_vertex, 'center'), 1)
        other_int = self._from['center'][from_vertex.next(self.sample_fol())]
        return self.get_oriented_edge(from_vertex, other_int, 'center',
                                      from_vertex.raw_endpoint(RIGHT, self.sample_fol()),
                                      self.sample_fol())


    def get_oriented_edge(self, from_vertex, to_vertex, edge_type,
                          crossing = None, fol = None):
        v = [from_vertex, to_vertex]
        candidates = [OrientedEdge((v[i], v[(i+1)%2], edge_type),
                                   (-1)**i)
                      for i in range(2)
                      if v[i] in self._from[edge_type] and
                      self._from[edge_type][v[i]] == v[(i+1)%2]]

        if len(candidates) == 1:
            return candidates[0]
        elif len(candidates) == 2:
            # this case may only happen if the foliation is two-sided
            assert(not self.sample_fol().is_bottom_side_moebius())
            bottom_interval = from_vertex if from_vertex.side == BOTTOM else to_vertex
            if (crossing <= bottom_interval.endpoint(RIGHT, fol)) == \
               (from_vertex.side == TOP):
                return candidates[0]
            return candidates[1]
            
        assert(False)

    def path_to_vector(self, path, is_signed):
        result = [0] * len(self._index_to_edge)
        # print self._edge_to_index
        # print path
        for oe in path:
            result[self._edge_to_index[oe.edge]] += oe.direction \
                                    if is_signed == SIGNED else 1
        return result

    def path_to_vector_teich(self, path, edge_index_to_elevation):
        edge_coeffs = [0] * len(self._index_to_edge)
        elevation = 1
        for oe in path:
            index = self._edge_to_index[oe.edge]
            if oe.direction == -1:
                elevation /= edge_index_to_elevation[index]
            edge_coeffs[index] += elevation
            if oe.direction == 1:
                elevation *= edge_index_to_elevation[index]

        return (vector(edge_coeffs), elevation)

    def _latex_(self):
        r"""
        Returns the LaTeX/TikZ representation of the train track.
        
        """
        from foliation_latex import FoliationLatex
        return FoliationLatex(self.sample_fol()).tikz_picture(
            train_tracks = [self])
        

    def get_symmetries_from(self, from_tt = None):
        self_map = (from_tt == None)
        if self_map:
            if hasattr(self, '_self_symmetries'):
                return self._self_symmetries
            from_tt = self

            
        if from_tt.invariants() != self.invariants():
            return []

        symmetries = []
        for interval in self.sample_fol().intervals():
            for hdir in [LEFT, RIGHT]:
                tt_map = self.transform(interval, hdir)
                if tt_map.domain == from_tt:
                    symmetries.append(tt_map)
        if self_map:
            self._self_symmetries = symmetries
        return symmetries


        
class TrainTrackPath(list):
    def reversed(self):
        return TrainTrackPath([oe.reversed() for 
                        oe in reversed(self)])

    @classmethod
    def from_edge(self, edge):
        return TrainTrackPath([OrientedEdge(edge, 1)])

    def start(self):
        return self[0].start()

    def end(self):
        return self[-1].end()


