r"""
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

from collections import namedtuple
from pseudo_anosov import PseudoAnosov
from foliation import SymmetryCoding
from transverse_curve import RestrictionCoding, get_codings, TransverseCurve
from base import RestrictionError, SaddleConnectionError

# ALLOWED_ERROR_OF_FOLIATIONS = 0.00001

PACandidate = namedtuple("PACandidate", "foliation, coding")

def find_pseudo_anosovs(foliation, depth):
    candidates = find_pseudo_anosov_candidates(foliation, depth)
    # print "BAD CANDIDATES: \n"
    # for x in candidates:
    #     if not is_really_pseudo_anosov(x):
    #         print x, '\n'
    result = []
    for cand in candidates:
        pa = create_pseudo_anosov(cand)
        if pa != None:
            result.append(pa)

    return result

# from train_track_map import tt_map_from_codings

def tt_map_from_codings_and_fol(start_fol, coding_list):
    from train_track_map import TrainTrackMap
    from foliation import SymmetryCoding
    from transverse_curve import RestrictionCoding, TransverseCurve
    tt_map_so_far = start_fol.train_track().identity_map()
    fol = start_fol
    for coding in coding_list:
        if isinstance(coding, SymmetryCoding):
            fol, tt_map = fol.transform(*coding)
        elif isinstance(coding, RestrictionCoding):
            tc = TransverseCurve(fol, coding.coding)
            fol, tt_map = tc.new_foliation()
        else:
            assert(False)

        tt_map_so_far = tt_map_so_far * tt_map            

    return tt_map_so_far

    
def create_pseudo_anosov(pa_cand):
    try:
        tt_map = tt_map_from_codings_and_fol(pa_cand.foliation, pa_cand.coding)
    except RestrictionError as ex:
        print "RestrictionError: ", ex
        return None
    except SaddleConnectionError as ex:
        print "SaddleConnectionError: ", ex
        return None

    try:
        return PseudoAnosov(tt_map)
    except ValueError:
        return None

    # if not tt_map.codomain.foliation.equals(tt_map.domain.foliation,
    #                            ALLOWED_ERROR_OF_FOLIATIONS):
    #     return None


            


def find_pseudo_anosov_candidates(foliation, depth,
                                  tt_map_so_far = None,
                                  coding_so_far = None):
    from train_track import TrainTrack
    from sage.rings.real_double import RDF
    # from foliation import SaddleConnectionError

    if tt_map_so_far == None:
        tt_map_so_far = foliation.train_track().identity_map()
        coding_so_far = []

    result = []

    for hdir in [RIGHT, LEFT]:
        for interval in foliation.intervals():
            final_fol, tt_map = foliation.transform(interval, hdir)
            final_tt_map = tt_map_so_far * tt_map
            final_coding = coding_so_far + [SymmetryCoding(interval,
                                                           hdir)]
            
            if final_tt_map.domain.sample_fol().permutation() ==\
               final_tt_map.codomain.sample_fol().permutation():
                # print final_tt_map.domain.foliation
                # print final_tt_map.codomain.foliation
                # print final_coding

                # BUG: this is buggy when the twist is longer than
                # the first interval
                # if final_tt_map.is_pseudo_anosov():
                lengths = final_tt_map.get_PF_interval_lengths(RDF)
                if lengths == None:
                    continue

                # m = final_tt_map.small_matrix().transpose()
                # try:
                #     eigenvalue, eigenvector = pf_eigen_data(m, RDF)
                # except NoPFEigenvectorError as ex:
                #     # print "No PF Eigendata: ", ex
                #     continue
                # if final_tt_map.is_pseudo_anosov():
                
                # print "Charpoly:", m.charpoly().factor()
                # print matrix(m, RDF).eigenvectors_right()
                # print "Stretch factor:", eigenvalue
                # print "Eigenvector:", eigenvector

                # f = final_tt_map.codomain.foliation
                try:
                    new_fol = final_fol.with_changed_lengths(lengths)
                except SaddleConnectionError:
                    # sometimes the generated new lengths lead to a foliation
                    # with a saddle connection
                    continue
                # ps = PseudoAnosov((tt_map, almost_final_tt_map, final_tt_map, m.charpoly(), eigenvalue, new_fol), final_coding)
                # ps = PseudoAnosov((m.charpoly(), eigenvalue, new_fol), final_coding)
                # ps = PseudoAnosov(new_fol, final_coding, m.charpoly(), eigenvalue,
                #                   final_tt_map)
                ps = PACandidate(new_fol, final_coding)
                result.append(ps)
                # result.extend([m, final_tt_map.edge_matrix().transpose()])



    if depth == 0:
        return result

    codings = get_codings(foliation)
    for c in codings:
        try:
            tc = TransverseCurve(foliation, c)
            new_fol, tt_map = tc.new_foliation()
            # m1 = tt_map.edge_matrix()
            # m2 = (tt_map_so_far * tt_map).edge_matrix()
            result.extend(find_pseudo_anosov_candidates(\
                        new_fol, depth - 1,
                        tt_map_so_far * tt_map,
                        coding_so_far + [RestrictionCoding(foliation, c)]))
                        # coding_so_far + [c, tt_map]))
                        # coding_so_far + [c, foliation, m1, m2]))
        except (RestrictionError, SaddleConnectionError):
            pass
    

    return result
    
def find_PAs_with_sing_type(sing_type, stop_after_time = 10, depth = 1):
    f = foliation_with_sing_type(sing_type)
    pas = []
    import time
    from base import SaddleConnectionError
    start = time.time()
    while time.time() - start < stop_after_time:
        f = f.with_changed_lengths()
        new_pas = find_pseudo_anosovs(f,depth)
        pas.extend(new_pas)
        # if len(new_pas) > 0:
        print str(len(new_pas)) + "pAs found."
    padict = {round(x.stretch_factor(), 6) : x for x in pas if
              x.stretch_factor() < 10}
    return padict
    
def foliation_with_sing_type(sing_type):
    assert(all(x%2 == 0 for x in sing_type))
    num_labels = sum(sing_type)/2
    labels = 2*range(num_labels)
    import random
    limit = 1000
    while limit > 0:
        limit -= 1
        random.shuffle(labels)
        f = Foliation(labels, 'moebius')
        if f.singularity_type() != sing_type:
            continue
        try:
            tt = f.train_track()
            return f
        except SaddleConnectionError:
            continue
    raise ValueError("No foliation with such singularity type has been"
                     "found.")
