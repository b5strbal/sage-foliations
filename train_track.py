from sage.structure.sage_object import SageObject
from mymath import mod_one
from collections import namedtuple
from sage.matrix.constructor import matrix, vector
from sage.rings.integer_ring import ZZ
from sage.symbolic.ring import var, SR
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

        self._coboundary_map = [None, None]
        self._tree_edge_indices = None

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

    def vertex_to_index(self, vertex):
        return self._vertex_to_index[vertex]


    def edges(self):
        return self._index_to_edge # list

    def kernel_from_singularities(self):
        circles = self.foliation.paths_around_singularities
        return [self.path_to_vector(circle, SIGNED) for circle in circles]

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

    def tree_edge_indices(self):
        if self._tree_edge_indices == None:
            m = self.coboundaries()
            self._tree_edge_indices = m.pivots()

        return self._tree_edge_indices

    def coboundaries(self):
        m = self.coboundary_map_from_vertices(COHOMOLOGY)
        m = m.transpose().echelon_form(include_zero_rows = False)
        return m
        

    def basis_of_cohomology(self, punctured):
        # rows generate the kernel
        if punctured:
            M = matrix.identity(ZZ, len(self.edges()))
        else:
            M = matrix(self.kernel_from_singularities()).right_kernel_matrix()
        U = self.vertex_map_normalizer()

        F = M * U.transpose()
        
        # rank of the vertex-to-edge matrix
        k = len(self.vertices()) - 1
        
        # The last k columns of F are killed by the image of c
        # after the normalization of U. 
        # So we just have to echelonize, and drop all the rows where non-zero
        # elements appear only in the last k columns, and renormalize.
        F.echelonize()
        # print F, '\n'
        M = F * U.transpose().inverse()

        # cutting off unnecessary rows
        col = row = 0
        while row < F.nrows():
            while col < F.ncols() - k and F[row,col]==0:
                col += 1
            if col == F.ncols() - k:
                break
            col += 1
            row += 1
        basis = matrix(M.rows()[:row])
        # coboundaries = matrix(M.rows()[row:])
        # coboundaries.echelonize()
        # print basis, '\n'
        # print coboundaries, '\n'
        # print self.coboundary_map_from_vertices('cohomology').transpose().echelon_form(), '\n'

        # What's left to do is to use the coboundaries to kill all non-zero entries in the
        # columns corresponding to the tree edge indices

        X = basis.matrix_from_columns(self.tree_edge_indices())
        to_subtract = X * self.coboundaries()
        
        return basis - to_subtract

    def vertex_map_normalizer(self):
        if not hasattr(self, '_U'):
            C = self.coboundary_map_from_vertices(COHOMOLOGY)
            D, U, V = C.smith_form() # D = U * C * V

            # Modifying U in the smith form such that D is a diagonal matrix
            # where the non-zero rows are in the end instead of at the beginning
            self._U = matrix(list(reversed(U.rows())))

        return self._U
        
    def small_vector_size(self):
        numints = len(self._index_to_vertex) / 2
        return numints if self.foliation.is_bottom_side_moebius() \
            else numints + 1
                    
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
        identity_size = X.nrows() if self._foliation.is_bottom_side_moebius() \
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
        return FoliationLatex(self._foliation).tikz_picture(
            train_tracks = [self])
        
        
    class Path(list):
        def reversed(self):
            return TrainTrack.Path([oe.reversed() for 
                            oe in reversed(self)])

        @classmethod
        def from_edge(self, edge):
            return TrainTrack.Path([OrientedEdge(edge, 1)])

        def start(self):
            return self[0].start()

        def end(self):
            return self[-1].end()



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
               
            new_map = TrainTrack.Map(domain = other.domain,
                                  codomain = self.codomain,
                                  vertex_map = new_vertex_map,
                                  edge_map = new_edge_map)
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
            return TrainTrack.Path(new_path)
            
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

        @classmethod
        def identity(cls, foliation):
            tt = foliation.train_track
            return TrainTrack.Map(tt, tt,
                                  {v:v for v in tt.vertices()},
                {e:TrainTrack.Path.from_edge(e) for e in tt.edges()})
                
        def is_self_map(self):
            return self.domain == self.codomain

        def action_on_cohomology(self, punctured = False):
            # when the surface is nonorientable, this matrix should be
            # shrinking with eigenvalue 1/alpha, since it is a pullback.
            basis = self.codomain.basis_of_cohomology(punctured)
            new_vectors = basis * self.edge_matrix(SIGNED).transpose()
            U = self.domain.vertex_map_normalizer()
            new_vectors *= U.transpose()
            # the rows represent the range of the linear map, so we transpose it
            return new_vectors.matrix_from_columns(range(new_vectors.nrows())).transpose()

        def invariant_cohomology(self, punctured = False):
            assert(self.is_self_map())
            action = self.action_on_cohomology(punctured)
            m = action - matrix.identity(ZZ, action.nrows())
            kernel = m.right_kernel_matrix()
            basis = self.domain.basis_of_cohomology(punctured)
            
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
                self.codomain.small_vector_size()))
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

        

        def teichmuller_polynomial(self, punctured):
            # from sage.rings.polynomial.polynomial_ring_constructor import \
            #     PolynomialRing
            from mymath import fixed_charpoly
            tt = self.domain
            inv_coh = self.invariant_cohomology(punctured)
            variables = get_variables(inv_coh.nrows())
            edge_index_to_elevation = [multi_power(variables, column)
                                       for column in inv_coh.columns()]
                                                   
            # ring = PolynomialRing(ZZ, variables)
            edge_matrix = matrix(SR, len(tt.edges()))
            vertex_matrix = matrix(SR, len(tt.vertices()))
            vertex_elevations = [None] * len(tt.vertices())

            # calculating the image of on vertex which uniquely defines
            # the lift of the train track map

            first_path = self.edge_map[tt.edges()[0]]
            vertex_matrix[0, tt.vertex_to_index(first_path.start())] = 1
            vertex_elevations[0] = 1
            # edge_matrix[0], elevation = tt.path_to_vector_teich(first_path,
            #                                                     edge_index_to_elevation)
            # print edge_index_to_elevation, '\n'
            while True:
                # find an edge such that the image of one endpoint is known
                i = find_next_edge_index(tt, vertex_matrix, edge_matrix)

                # if no edges are left, break out
                if i == -1:
                    break

                edge = tt.edges()[i]
                path = self.edge_map[edge]
                v, elevation = tt.path_to_vector_teich(path,
                                                       edge_index_to_elevation)
                # print 'path: ', path
                # print 'path elevation: ', elevation
                # print 'edge: ', edge
                # print 'edge elevation: ', edge_index_to_elevation[i]
                vindex = [tt.vertex_to_index(edge[j]) for j in [0,1]]
                # print elevation
                net_elev = elevation / edge_index_to_elevation[i]
                # print net_elev, (-1)**0, (-1)**1, net_elev ** (-1), net_elev**1
                endpoint_ind = [path.start(), path.end()]
                endpoint_ind = [tt.vertex_to_index(x) for x in endpoint_ind]
                for j in [0,1]:
                    if vertex_elevations[vindex[j]] == None:
                        mult_by = net_elev if j == 1 else 1/net_elev
                        vertex_elevations[vindex[j]] = \
                                vertex_elevations[vindex[(j+1)%2]] * mult_by
                        vertex_matrix[vindex[j], endpoint_ind[j]] = \
                                            vertex_elevations[vindex[j]]
                        break
                    
                v *= vertex_elevations[vindex[0]]
                # print v
                # print edge_matrix
                edge_matrix[i] = v
                # print edge_matrix, '\n'
                # print vertex_matrix, '\n'

            # print edge_matrix, '\n'
            # print vertex_matrix, '\n'
            # return (edge_matrix, vertex_matrix)
            p1 = fixed_charpoly(edge_matrix, variables)
            p2 = fixed_charpoly(vertex_matrix, variables)
            # print p1, p2
            # return (p1, p2)
            # print SR(p1), SR(p2)

            tpoly = (SR(p1)/SR(p2)).simplify_rational()
            # tpoly.reduce()
            return tpoly
                    
                
        def is_pseudo_anosov(self):
            return is_perron_frobenius(self.edge_matrix())



def find_next_edge_index(tt, vertex_matrix, edge_matrix):
    for i in xrange(edge_matrix.nrows()):
        if not edge_matrix[i].is_zero():
            continue
        edge = tt.edges()[i]
        for k in [0,1]:
            vindex = tt.vertex_to_index(edge[k])
            if not vertex_matrix[vindex].is_zero():
                return i
    return -1
    


def reduce(tt):
    """Create a projection to a traintrack without 2-prong vertices."""
    
    pass


def get_variables(n):
    if n == 0:
        return []
    names = ['u','t','s','v','w','y','z']
    if n == 1:
        return [var('u')]
    if n <= 7:
        return var(names[:n])
    names = ['x' + str(i) for i in xrange(n)]
    return var(names)

def multi_power(variables, powers):
    from sage.misc.misc_c import prod
    assert(len(powers) == len(variables))
    if len(powers) == 0:
        return 1
    return prod(variables[i]**powers[i] for i in xrange(len(powers)))
