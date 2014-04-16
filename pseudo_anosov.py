from sage.structure.sage_object import SageObject
from sage.symbolic.ring import var, SR
from sage.rings.real_double import RDF
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix
from mymath import fixed_charpoly, NoPFEigenvectorError, pf_eigen_data
from constants import epsilon, TRAIN_TRACK_MODULE
from interval import Interval

class PseudoAnosov(SageObject):
    def __init__(self, tt_map):
        if not tt_map.is_self_map():
            raise ValueError("The train track map should be a self-map.")
        m = tt_map.edge_matrix().transpose()
        try:
            eigenvalue, eigenvector = pf_eigen_data(m, RDF)
        except NoPFEigenvectorError as ex:
            # print "No PF Eigendata: ", ex
            raise ValueError("The transition map of the train track map "
                             "doesn't have a PF eigenvector.")

        tt = tt_map.domain

        # making sure that the vertex conditions are met
        C = tt.coboundary_map_from_vertices(TRAIN_TRACK_MODULE).transpose()
        if abs(C*eigenvector) > epsilon:
            raise ValueError("The PF eigenvector is not an allowable "
                             "weights vector.")

        self._tt_map = tt_map
        self._tt_weights = eigenvector
        self._stretch_factor = eigenvalue

        k = len(eigenvector)/3
        lengths = list(eigenvector[:k])
        if not tt.sample_fol().is_bottom_side_moebius():
            edge = tt.get_edge_from(Interval(1,0), 'center')
            n = edge[1].index
            twist = sum(eigenvector[k:k+1+n])
            lengths.append(twist)
        
        # this might raise a SaddleConnectionError, but it should be very rare
        # as the map is pseudo-anosov
        self._stable_fol = tt.sample_fol().with_changed_lengths(lengths)

    def _repr_(self):
        s = "Pseudo-anosov map with stretch factor "
        s += str(self.stretch_factor())
        s += " and characteristic polynomial "
        s += str(self.teichmuller_polynomial(False))
        return s

    def stretch_factor(self):
        return self._stretch_factor

    


    def teichmuller_polynomial(self, is_punctured):
        # from sage.rings.polynomial.polynomial_ring_constructor import \
        #     PolynomialRing
        tt = self._tt_map.domain
        tt_map = self._tt_map
        inv_coh = self._tt_map.invariant_cohomology(is_punctured)
        variables = get_variables(inv_coh.nrows())
        edge_index_to_elevation = [multi_power(variables, column)
                                   for column in inv_coh.columns()]

        # ring = PolynomialRing(ZZ, variables)
        edge_matrix = matrix(SR, len(tt.edges()))
        vertex_matrix = matrix(SR, len(tt.vertices()))
        vertex_elevations = [None] * len(tt.vertices())

        # calculating the image of on vertex which uniquely defines
        # the lift of the train track map

        first_path = tt_map.edge_map[tt.edges()[0]]
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
            path = tt_map.edge_map[edge]
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
    assert(len(powers) == len(variables))
    if len(powers) == 0:
        return 1
    return prod(variables[i]**powers[i] for i in xrange(len(powers)))

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
    
