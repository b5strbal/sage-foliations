from separatrix import Separatrix
from arc import Arc
from collections import namedtuple
from sage.structure.sage_object import SageObject

class RestrictionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


Coding = namedtuple("Coding", "side, index, end, num_flips1, num_flips2")


class TransverseCurve(SageObject):
    def __init__(self, foliation, coding):
        interval = foliation.interval(coding.side, coding.index)

        # Checking if the curve is transverse
        if (coding.num_flips1 % 2 == 1) != \
           (coding.num_flips2 % 2 == 1) != \
                    interval.is_flipped():
            raise RestrictionError("Curve is not transverse to the foliation")


        other_int = interval.pair()
        other_end = coding.end
        if interval.is_flipped():
            other_end = (other_end + 1) % 2
        sep1 = Separatrix(foliation, interval, end = coding.end,
                          number_of_flips_to_stop = coding.num_flips1)
        sep2 = Separatrix(foliation, other_int, end = other_end,
                          number_of_flips_to_stop = coding.num_flips2)

        if sep1.is_flipped() == (coding.end == 1):
            openness = ('open','closed')
        else:
            openness = ('closed','open')

        self._sep = [sep1, sep2]
        arcs = [Arc(self._sep[i].endpoint, self._sep[(i+1)%2].endpoint,
                    *openness) for i in range(2)]
            

        self._arc = Arc(sep1.endpoint, sep2.endpoint, *openness)
        import itertools
        intersections = itertools.chain(sep1.intersections()[:-1],
                                        sep2.intersections()[:-1])

        # if one of the separatrices is longer, both arcs are considered,
        # so if the current one is wrong, we try the other one
        if coding.num_flips1 > 0 or coding.num_flips2 > 0:
            x = next(intersections)
            if self._arc.contains(x):
                self._arc = Arc(sep2.endpoint, sep1.endpoint, *openness)
                self._sep = [sep2, sep1]

        if self._arc.length() > 0.5 and sep1.end_side == sep2.end_side:
            raise RestrictionError("The curve is one-side and has length"
                                   " greater than half, so we can't compute"
                                   " the train track transition map. ")

        for x in intersections:
            if self._arc.contains(x):
                raise RestrictionError("The curve is self-intersecting")

        self._direction = 'right' if openness[0] == 'closed' else 'left'
        self._coding = coding

    def _repr_(self):
        s = ''
        s += 'Arc: ' + repr(self._arc) + '(' + repr(self._arc.length()) + ')'
        return s

    def coding(self):
        return self._coding
                           
    def is_one_sided(self):
        return self._sep[0].end_side == self._sep[1].end_side

    def new_foliation(self):
        foliation = self._sep[0].foliation
        separatrices = Separatrix.get_all(foliation, self._arc)
        from transition_map import new_foliation
        starting_side = self._sep[0].end_side if self._direction == 'right' \
                       else self._sep[1].end_side
        return new_foliation(separatrices, self._sep[0].endpoint,
                             starting_side,
                             is_one_sided = self.is_one_sided(),
                             direction = self._direction,
                             ending_point = self._sep[1].endpoint)



def get_codings(foliation):
    codings = []
    for interval in foliation.intervals():
        if not interval.is_flipped():
            codings.append(Coding(interval.side, interval.index,
                                  0, 0, 0))
        else:
            codings.append(Coding(interval.side, interval.index,
                                  0, 0, 1))
            codings.append(Coding(interval.side, interval.index,
                                  1, 0, 1))
                                
    return codings

PseudoAnosov = namedtuple("PseudoAnosov", "foliation, coding") 

def find_pseudo_anosovs(foliation, depth):
    candidates = find_pseudo_anosov_candidates(foliation, depth)
    return candidates

def find_pseudo_anosov_candidates(foliation, depth,
                                  tt_map_so_far = None,
                                  coding_so_far = None):
    from train_track import TrainTrack
    from mymath import NoPFEigenvectorError, pf_eigen_data
    if tt_map_so_far == None:
        tt_map_so_far = TrainTrack.Map.identity(foliation)
        coding_so_far = []

    result = []

    for is_reversed in [False, True]:
        if is_reversed:
            almost_final_fol, tt_map = foliation.reversed()
            almost_final_tt_map = tt_map_so_far * tt_map
            almost_final_coding = coding_so_far + [is_reversed]
        else:
            almost_final_fol = foliation
            almost_final_coding = coding_so_far
            almost_final_tt_map = tt_map_so_far
        for rotateby in range(foliation.num_intervals(0)):
            final_fol, tt_map = almost_final_fol.rotated(rotateby)
            final_tt_map = almost_final_tt_map * tt_map
            final_coding = almost_final_coding + [rotateby]
            # print final_coding
            # print final_tt_map.domain.foliation
            # print final_tt_map.codomain.foliation
            if final_tt_map.domain.foliation.permutation() ==\
               final_tt_map.codomain.foliation.permutation():
                # print final_tt_map.domain.foliation
                # print final_tt_map.codomain.foliation
                # print final_coding
                m = final_tt_map.small_matrix().transpose()
                try:
                    eigenvalue, eigenvector = pf_eigen_data(m, RDF)
                except NoPFEigenvectorError as ex:
                    print "No PF Eigendata: ", ex
                    continue
                # if final_tt_map.is_pseudo_anosov():
                
                # print "Charpoly:", m.charpoly().factor()
                # print matrix(m, RDF).eigenvectors_right()
                # print "Stretch factor:", eigenvalue
                # print "Eigenvector:", eigenvector

                f = final_tt_map.codomain.foliation
                try:
                    new_fol = f.with_changed_lengths(list(eigenvector))
                except SaddleConnectionError:
                    # sometimes the generated new lengths lead to a foliation
                    # with a saddle connection
                    continue
                ps = PseudoAnosov((m.charpoly(), eigenvalue, new_fol), final_coding)
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
                        coding_so_far + [c]))
                        # coding_so_far + [c, foliation, m1, m2]))
        except RestrictionError: # also SaddleConnectionError?
            pass
    

    return result
    
