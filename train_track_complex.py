from sage.structure.sage_object import SageObject
from train_track_class import TrainTrackClass
from transverse_curve import Coding, TransverseCurve

class TrainTrackComplex(DiGraph):
    def __init__(self, first_tt):
        # adding the first vertex
        super(TrainTrackComplex, self).__init__(multiedges = True,
                                                loops = True)
        self.get_present_vertex(first_tt)
        self._total_tt_map_searched = 0

    # def add_vertex(self, tt):
    #     tt_class = TrainTrackClass(tt)
    #     self._vertices.add(tt_class)

    def add_tt_map(self, tt_map):
        self._total_tt_map_searched += 1
        start = self.get_present_vertex(tt_map.codomain)
        end = self.get_present_vertex(tt_map.domain)

        # print tt_map.codomain
        # print start
        # print tt_map.codomain.get_symmetries(start.tt_repr())
        # print tt_map.domain
        # print end
        # print end.tt_repr().get_symmetries(tt_map.domain)

        # print tt_map
        # print tt_map.codomain.get_symmetries(start.tt_repr())[0], '\n'
        # print tt_map.codomain.get_symmetries(start.tt_repr())[0] *\
        #              tt_map, '\n'
        # print tt_map *\
        #              end.tt_repr().get_symmetries(tt_map.domain)[0]

        new_tt_map = start.tt_repr().get_symmetries_from(tt_map.codomain)[0]*\
                     tt_map *\
                     tt_map.domain.get_symmetries_from(end.tt_repr())[0]

        print new_tt_map
        for edge in self.edges_incident():
            if edge[2] == new_tt_map:
                return

        print 'ufff'
        self.add_edge((start, end, new_tt_map))
        print (start, end, new_tt_map)
        
    def get_present_vertex(self, tt):
        tt_class = TrainTrackClass(tt)
        for v in self.vertices():
            if v == tt_class:
                return v
        self.add_vertex(tt_class)
        return tt_class


    def build(self, stop_after_time = 10):
        import random, time
        from myexceptions import SaddleConnectionError
        start = time.time()
        total_count_at_start = len(self.edges())
        count = 0
        while time.time() - start < stop_after_time:
            v = random.choice(self.vertices())

            # new foliation with random lengths
            try:
                fol = v.tt_repr().sample_fol().with_changed_lengths()
            except SaddleConnectionError:
                continue

            codings = get_codings(fol)
            for c in codings:
                try:
                    tc = TransverseCurve(fol, c)
                except RestrictionError:
                    continue

                try:
                    new_fol, tt_map = tc.new_foliation()
                except SaddleConnectionError:
                    continue

                self.add_tt_map(tt_map)
                count += 1

        print "tt_maps considered in this call: " + repr(count)
        print "Of these, " + repr(len(self.edges()) - total_count_at_start) + " were new."
        print "tt_maps considered overall: " + repr(self._total_tt_map_searched)
        print "Of these, " + repr(len(self.edges())) + " different ones."

    def search_from_vertex(self, depth, current_vertex, tt_map_so_far = None):
        from pseudo_anosov import PseudoAnosov
        if tt_map_so_far == None:
            tt_map_so_far = current_vertex.tt_repr().identity_map()

        if depth == 0:
            try:
                return [PseudoAnosov(tt_map_so_far)]
            except ValueError as ex:
                print ex
                return []

        result = []
        for edge in self.edges_incident(current_vertex):
            result.extend(self.search_from_vertex(depth - 1,
                                                  edge[1],
                                                  tt_map_so_far*edge[2]))
        return result

    def search(self, depth):
        result = []
        for d in range(1, depth + 1):
            for v in self.vertices():
                result.extend(self.search_from_vertex(d, v))
        return result


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




from examples import family_A_foliation
f = family_A_foliation(3, 1, False)
tt = f.train_track
ttc = TrainTrackComplex(tt)















# class TrainTrackComplex(SageObject):
#     def __init__(self):
#         self._vertices = set()
#         self._edges = {}

#     def add_vertex(self, tt):
#         tt_class = TrainTrackClass(tt)
#         if tt_class not in self._vertices:
#             self._vertices.add(tt_class)
#             self._edges[tt_class] = set()

#     def add_edge(self, tt_map):
#         self.add_vertex(tt_map.codomain)
#         self.add_vertex(tt_map.domain)
#         start = TrainTrackClass(tt_map.codomain)
#         self._edges[start].add(tt_map)
    
    
        
