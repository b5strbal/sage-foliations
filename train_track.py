from sage.structure.sage_object import SageObject
from sage.graphs.digraph import DiGraph
from foliation import mod_one
from collections import namedtuple

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
        done = set()
        for interval in foliation.intervals():
            if not interval in done:
                pair = interval.pair()
                self._from['pair'][interval] = pair
                done.add(pair)

        for interval in foliation.intervals():
            if not foliation.is_bottom_side_moebius():
                otherside = (interval.side + 1) % 2
                containing_int = foliation.in_which_interval(\
                        interval.endpoint(0), otherside)
            else:
                containing_int = foliation.in_which_interval(\
                        mod_one(interval.endpoint(0) + 0.5), 0)
            self._from['center'][interval] = containing_int


    @property
    def foliation():
        return self._foliation

    def get_edge_from(self, from_vertex, edge_type):
        return (from_vertex, self._from[edge_type][from_vertex], edge_type)

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
            top_interval = from_vertex if from_vertex.side == 0 else to_vertex
            if (crossing < top_interval.midpoint()) == (from_vertex.side == 0):
                return candidates[0]
            return candidates[1]
            
        assert(False)
            

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
            new_vertex_map = {v:self.vertex_map[other.vertex_map[v]]
                              for v in other.vertex_map}
            new_edge_map = {e:self._map_path[other.edge_map[e]]
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
            
        
                                  
                                
                                  


