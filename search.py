from collections import namedtuple
from transverse_curve import Coding, TransverseCurve

def get_codings(foliation):
    codings = []
    for interval in foliation.intervals():
        if not interval.is_flipped(foliation):
            codings.append(Coding(interval.side, interval.index,
                                  0, 0, 0))
        else:
            codings.append(Coding(interval.side, interval.index,
                                  0, 0, 1))
            codings.append(Coding(interval.side, interval.index,
                                  1, 0, 1))
                                
    return codings

# PseudoAnosov = namedtuple("PseudoAnosov", "foliation, coding, charpoly, eigenvalue,"
#                           "tt_map") 
# PseudoAnosov = namedtuple("PseudoAnosov", "tt_map, coding")

# PACandidate = namedtuple("PACandidate", "foliation, coding")

# def find_pseudo_anosovs(foliation, depth):
#     candidates = find_pseudo_anosov_candidates(foliation, depth)
#     # print "BAD CANDIDATES: \n"
#     # for x in candidates:
#     #     if not is_really_pseudo_anosov(x):
#     #         print x, '\n'
#     result = []
#     for cand in candidates:
#         pa = create_pseudo_anosov(cand)
#         if pa != None:
#             result.append(pa)

#     return result

# ALLOWED_ERROR_OF_FOLIATIONS = 0.00001

# def print_charpoly(tt_map):
#     m = tt_map.small_matrix()
#     if m.is_square():
#         print m.charpoly()
#     else:
#         print "Small matrix not square."
    
def tt_map_from_codings(foliation, coding_list):
    from train_track_map import TrainTrackMap
    tt_map_so_far = TrainTrackMap.identity(foliation.train_track)
    for coding in coding_list:
        # the last two elements of the list are the reversing and the rotation
        if not isinstance(coding, Coding):
            break
        tc = TransverseCurve(foliation, coding)
        foliation, tt_map = tc.new_foliation()
        tt_map_so_far = tt_map_so_far * tt_map
        # print_charpoly(tt_map_so_far)
        # print tt_map_so_far.domain.basis_of_cohomology(False), '\n'
        
    if coding_list[-2]: # reversed
        foliation, tt_map = foliation.reversed()
        tt_map_so_far = tt_map_so_far * tt_map
        # print_charpoly(tt_map_so_far)
        # print tt_map_so_far.domain.basis_of_cohomology(False), '\n'
    
        

    foliation, tt_map = foliation.rotated(coding_list[-1])
    tt_map_so_far = tt_map_so_far * tt_map
    # print tt_map_so_far.domain.basis_of_cohomology(False), '\n'
    # print_charpoly(tt_map_so_far)
    # print '\n'
    return tt_map_so_far

# def create_pseudo_anosov(pa_cand):
#     try:
#         tt_map = tt_map_from_codings(pa_cand.foliation, pa_cand.coding)
#     except RestrictionError as ex:
#         print "RestrictionError: ", ex
#         return None
#     except SaddleConnectionError as ex:
#         print "SaddleConnectionError: ", ex
#         return None
        
#     if not tt_map.codomain.foliation.equals(tt_map.domain.foliation,
#                                ALLOWED_ERROR_OF_FOLIATIONS):
#         return None

#     return PseudoAnosov(tt_map, pa_cand.coding)
            

def find_pseudo_anosovs(foliation, depth,
                                  tt_map_so_far = None,
                                  coding_so_far = None):
    from train_track_map import TrainTrackMap
    # from mymath import NoPFEigenvectorError, pf_eigen_data
    # from sage.rings.real_double import RDF
    from myexceptions import SaddleConnectionError

    if tt_map_so_far == None:
        tt_map_so_far = TrainTrackMap.identity(foliation.train_track)
        coding_so_far = []

    result = []

    for is_reversed in [False, True]:
        if is_reversed:
            almost_final_fol, tt_map = foliation.reversed()
            almost_final_tt_map = tt_map_so_far * tt_map
            almost_final_coding = coding_so_far + [is_reversed]
        else:
            almost_final_fol = foliation
            almost_final_coding = coding_so_far + [is_reversed]
            almost_final_tt_map = tt_map_so_far
        for rotateby in range(foliation.num_intervals(0)):
            final_fol, tt_map = almost_final_fol.rotated(rotateby)
            final_tt_map = almost_final_tt_map * tt_map
            final_coding = almost_final_coding + [rotateby]

            try:
                result.append(PseudoAnosov(final_tt_map))
            except (SaddleConnectionError, ValueError) as ex:
                print ex
                continue

            # # print final_coding
            # # print final_tt_map.domain.foliation
            # # print final_tt_map.codomain.foliation
            # if final_tt_map.domain.foliation.permutation() ==\
            #    final_tt_map.codomain.foliation.permutation():
            #     # print final_tt_map.domain.foliation
            #     # print final_tt_map.codomain.foliation
            #     # print final_coding

            #     # BUG: this is buggy when the twist is longer than
            #     # the first interval
            #     m = final_tt_map.small_matrix().transpose()
            #     try:
            #         eigenvalue, eigenvector = pf_eigen_data(m, RDF)
            #     except NoPFEigenvectorError as ex:
            #         # print "No PF Eigendata: ", ex
            #         continue
            #     # if final_tt_map.is_pseudo_anosov():
                
            #     # print "Charpoly:", m.charpoly().factor()
            #     # print matrix(m, RDF).eigenvectors_right()
            #     # print "Stretch factor:", eigenvalue
            #     # print "Eigenvector:", eigenvector

            #     f = final_tt_map.codomain.foliation
            #     try:
            #         new_fol = f.with_changed_lengths(list(eigenvector))
            #     except SaddleConnectionError:
            #         # sometimes the generated new lengths lead to a foliation
            #         # with a saddle connection
            #         continue
            #     # ps = PseudoAnosov((tt_map, almost_final_tt_map, final_tt_map, m.charpoly(), eigenvalue, new_fol), final_coding)
            #     # ps = PseudoAnosov((m.charpoly(), eigenvalue, new_fol), final_coding)
            #     # ps = PseudoAnosov(new_fol, final_coding, m.charpoly(), eigenvalue,
            #     #                   final_tt_map)
            #     ps = PACandidate(new_fol, final_coding)
            #     result.append(ps)
            #     # result.extend([m, final_tt_map.edge_matrix().transpose()])


    if depth == 0:
        return result

    codings = get_codings(foliation)
    for c in codings:
        try:
            tc = TransverseCurve(foliation, c)
        except RestrictionError:
            pass
            
        try:
            new_fol, tt_map = tc.new_foliation()
        except SaddleConnectionError:
            pass

            # m1 = tt_map.edge_matrix()
            # m2 = (tt_map_so_far * tt_map).edge_matrix()
        result.extend(find_pseudo_anosovs(\
                                          new_fol, depth - 1,
                                          tt_map_so_far * tt_map,
                                          coding_so_far + [c]))
        # coding_so_far + [c, tt_map]))
        # coding_so_far + [c, foliation, m1, m2]))
    

    return result
    
