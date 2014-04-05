from sage.structure.sage_object import SageObject
from mymath import mod_one
from collections import namedtuple
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from mymath import is_perron_frobenius
from constants import *

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
        self._foliation = foliation

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
                pair = interval.pair()
                self._from['pair'][interval] = pair
                edge = (interval, pair, 'pair')
                self._edge_to_index[edge] = edge_count
                self._index_to_edge.append(edge)
                edge_count += 1

        for interval in foliation.intervals():
            side = interval.endpoint_side(LEFT)
            containing_int = foliation.in_which_interval(\
                            interval.endpoint(LEFT), (side + 1) % 2)
            self._from['center'][interval] = containing_int
            edge = (interval, containing_int, 'center')
            self._edge_to_index[edge] = edge_count
            self._index_to_edge.append(edge)
            edge_count += 1

    def __eq__(self, other):
        # print self._index_to_edge
        # print other._index_to_edge
        # print self._index_to_edge == other._index_to_edge
        # print '\n'
        return self._index_to_edge == other._index_to_edge

    @property
    def foliation(self):
        return self._foliation

    def vertices(self):
        return self._index_to_vertex

    def edges(self):
        return self._index_to_edge # list

    def kernel_from_singularities(self):
        circles = self.foliation.paths_around_singularities
        return [self.path_to_vector(circle, signed = True) for circle in circles]

    def coboundary_map_from_vertices(self, cochain_type):
        m = matrix(ZZ, len(self._index_to_edge), len(self._index_to_vertex))
        if cochain_type == 'cohomology':
            for i in range(len(self._index_to_edge)):
                e = self._index_to_edge[i]
                m[i, self._vertex_to_index[e[0]]] -= 1
                m[i, self._vertex_to_index[e[1]]] += 1
        elif cochain_type == 'train track module':
            for i in range(len(self._index_to_edge)):
                e = self._index_to_edge[i]
                if e[2] == 'center':
                    x = 1
                else:
                    x = -1
                m[i, self._vertex_to_index[e[0]]] += x
                m[i, self._vertex_to_index[e[1]]] += x
            
        return m
         
    def small_vector_size(self):
        numints = len(self._index_to_vertex) / 2
        return numints if self.foliation.is_bottom_side_moebius() \
            else numints + 1
                    
    def matrix_to_reduce_dimension(self):
        if hasattr(self, '_reducing_matrix'):
            return self._reducing_matrix
        from copy import copy
        from sage.rings.rational import Rational
        X = self.coboundary_map_from_vertices('train track module')
        # flipping the matrix over vertically (so the echelonizing works
        # in the way we want) and transpose
        X = matrix(list(reversed(X.rows()))).transpose()
        X = X.echelon_form(include_zero_rows = False)
        
        # if the train track is non-orientable, the last row will be a
        # multiple of two, so we have to normalize it and re-echelonize.
        n = X.nrows()
        X = copy(X).with_row_set_to_multiple_of_row(n-1, n-1, Rational('1/2'))
        X = X.echelon_form()

        # now we need to cut off the identity part of the matrix and
        # append a small block of identity in the different direction
        X = matrix(list(reversed(X.rows())))
        X = -matrix(list(reversed(X.columns()[X.nrows():]))).transpose()
        X = matrix(matrix.identity(X.ncols()).rows() + X.rows())
        self._reducing_matrix = X
        return X

    def get_edge_from(self, from_vertex, edge_type):
        return (from_vertex, self._from[edge_type][from_vertex], edge_type)
        
    def get_center_edge(self, from_vertex, end):
        if end == 0:
            return OrientedEdge(self.get_edge_from(from_vertex, 'center'), 1)
        other_int = self._from['center'][from_vertex.next()]
        return self.get_oriented_edge(from_vertex, other_int, 'center',
                                      from_vertex.raw_endpoint(RIGHT))


    def get_oriented_edge(self, from_vertex, to_vertex, edge_type,
                          crossing = None):
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
            assert(not self._foliation.is_bottom_side_moebius())
            bottom_interval = from_vertex if from_vertex.side == BOTTOM else to_vertex
            if (crossing <= bottom_interval.endpoint(RIGHT)) == \
               (from_vertex.side == TOP):
                return candidates[0]
            return candidates[1]
            
        assert(False)

    def path_to_vector(self, path, signed):
        result = [0] * len(self._index_to_edge)
        # print self._edge_to_index
        # print path
        for oe in path:
            result[self._edge_to_index[oe.edge]] += oe.direction \
                                    if signed else 1
        return result


    def _latex_(self):
        r"""
        Returns the LaTeX/TikZ representaion of the train track.
        
        """
        from foliation_latex import FoliationLatex
        return FoliationLatex(self._foliation).tikz_picture(
            train_tracks = [self])
        
        
    class Path(list):
        def reversed(self):
            return TrainTrack.Path([oe.reversed() for 
                            oe in reversed(self)])

        @classmethod
        def from_edge(self, edge):
            return TrainTrack.Path([OrientedEdge(edge, 1)])

        # def __getitem__(self, index):
        #     result = list.__getitem__(self, index)
        #     try:
        #         return TrainTrack.Path(result)
        #     except TypeError:
        #         return result

        # def __add__(self, other):
        #     return TrainTrack.Path(self.oriented_edge_list +
        #                            other.oriented_edge_list)

               
    class Map(namedtuple("TrainTrackMap", "domain, codomain,"
                         "vertex_map, edge_map")):
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
               
            return TrainTrack.Map(domain = other.domain,
                                  codomain = self.codomain,
                                  vertex_map = new_vertex_map,
                                  edge_map = new_edge_map)
               
        def _map_path(self, path):
            new_path = []
            for oe in path:
                p = self.edge_map[oe.edge]
                if oe.direction == 1:
                    new_path.extend(p)
                else:
                    new_path.extend(p.reversed())
            return TrainTrack.Path(new_path)
            
        def _cohomology_kernel(self):
            # this only works if domain and codomain are combinatorically same
            # the underlying foliations need not be the same
            tt = self.domain
            
            m = matrix([tt.path_to_vector(self.edge_map[edge], signed = True)
                        for edge in tt.edges()])
            m = m.transpose() - matrix.identity(ZZ, m.nrows())
            mlist = m.rows()
            mlist.extend(tt.kernel_from_singularities())
            return matrix(mlist).right_kernel_matrix()

        @classmethod
        def identity(cls, foliation):
            tt = foliation.train_track
            return TrainTrack.Map(tt, tt,
                                  {v:v for v in tt.vertices()},
                {e:TrainTrack.Path.from_edge(e) for e in tt.edges()})
                
        def is_self_map(self):
            return self.domain == self.codomain

        def invariant_cohomology(self):
            # rows generate the kernel
            M = self._cohomology_kernel()
            C = self.domain.coboundary_map_from_vertices('cohomology')
            D, U, V = C.smith_form() # D = U * C * V

            # Modifying U in the smith form such that D is a diagonal matrix
            # where the non-zero rows are in the end instead of at the beginning
            U = matrix(list(reversed(U.rows())))

            F = M * U.transpose()
            k = C.rank() 
            # The first k columns of F are killed by the image of c
            # after the normalization of U. 
            # So we just have to echelonize, and drop all the rows where non-zero
            # elements appear only in the last k columns, and renormalize.
            F.echelonize()
            col = 0
            for row in range(F.nrows()):
                while col < F.ncols() - k and F[row,col]==0:
                    col += 1
                if col == F.ncols() - k:
                    break
                col += 1
            # cutting off unnecessary rows
            F = matrix(F.rows()[:row])
            if F.nrows() == 0:
                return F
            M = F * U.transpose().inverse()
            return M

        def small_matrix(self):
            # The domain and the codomain should be both one-sided or both
            # two-sided, otherwise the matrix won't be a square matrix
            # Usually this is not a problem, since we only call this method
            # is the underlying permutations are the same which is a much stronger
            # condition.
            m = self.domain.matrix_to_reduce_dimension()
            # print m
            result = m.transpose() * self.edge_matrix()
            result = result.matrix_from_columns(range(
                self.codomain.small_vector_size()))
            return result
                                  
        def edge_matrix(self):
            if not hasattr(self, '_edge_matrix'):
                tt = self.domain
                self._edge_matrix = matrix([self.codomain.path_to_vector(
                    self.edge_map[edge],
                    signed = False)
                                            for edge in self.domain.edges()])
            return self._edge_matrix

            
        def is_pseudo_anosov(self):
            return is_perron_frobenius(self.edge_matrix())



def reduce(tt):
    """Create a projection to a traintrack without 2-prong vertices."""
    
    pass
