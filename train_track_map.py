from constants import *
# from mymath import is_perron_frobenius
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



    # def small_matrix(self):
    #     # The domain and the codomain should be both one-sided or both
    #     # two-sided, otherwise the matrix won't be a square matrix
    #     # Usually this is not a problem, since we only call this method
    #     # is the underlying permutations are the same which is a much stronger
    #     # condition.
    #     m = self.domain.matrix_to_reduce_dimension()
    #     # print m
    #     # print self.domain.foliation
    #     # print self._edge_matrix, '\n'
    #     result = m.transpose() * self.edge_matrix()
    #     # print result, '\n'
    #     result = result.matrix_from_columns(range(
    #         self.codomain.small_vector_size()))
    #     # print result, '\n'

    #     return result

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




    # def is_pseudo_anosov(self):
    #     return is_perron_frobenius(self.edge_matrix())






def tt_map_from_codings(start_tt,coding_list):
    from train_track_map import TrainTrackMap
    from foliation import SymmetryCoding
    from transverse_curve import RestrictionCoding, TransverseCurve
    tt_map_so_far = start_tt.identity_map()
    for coding in coding_list:
        if isinstance(coding, SymmetryCoding):
            tt_map = tt_map_so_far.domain.transform(*coding)
        elif isinstance(coding, RestrictionCoding):
            tc = TransverseCurve(coding.foliation, coding.coding)
            tt_map = tc.new_foliation()[1]
        else:
            assert(False)

        tt_map_so_far = tt_map_so_far * tt_map            

    return tt_map_so_far


