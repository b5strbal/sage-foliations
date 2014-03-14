
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


class TrainTrack(DiGraph):
    def __init__(self, foliation):
        self._foliation = foliation
        DiGraph.__init__(self, multiedges = True)

        done = set()
        for interval in foliation.intervals():
            if not interval in done:
                pair = interval.pair()
                self.add_edge(interval, pair, 'pair')
                done.add(interval)
                done.add(pair)

        for interval in foliation.intervals():
            if not foliation.is_bottom_side_moebius():
                otherside = (interval.side + 1) % 2
                containing_int = foliation.in_which_interval(\
                        interval.endpoint(0), otherside)
            else:
                containing_int = foliation.in_which_interval(\
                        mod_one(interval.endpoint(0) + 0.5), 0)
            self.add_edge(interval, containing_int, 'center')

    @property
    def foliation():
        return self._foliation

    

    def get_oriented_edge(self, start, end, label, crossing = None):
        has_forward = self.has_edge(start, end, label)
        has_backwards = self.has_edge(end, start, label)
        assert(has_forward or has_backwards)
        if not has_forward:
            return OrientedEdge((end, start, label), -1)
        if not has_backwards:
            return OrientedEdge((start, end, label), 1)
        top_interval = start if start.side == 0 else end
        if (crossing < top_interval.midpoint()) == (start.side == 0):
            return OrientedEdge((start, end, label), 1)
        return OrientedEdge((end, start, label), -1)

    # def get_first_edge(self, separatrix, left_or_right):
    #     first_interval = separatrix.traversed_intervals[0]
    #     if left_or_right == 'left':
    #         first_interval = first_interval.prev()
    #     if self._foliation.is_bottom_side_moebius():
    #         return self.get_oriented_edge(first_interval, 
    #                 separatrix.traversed_intervals[3], 'center')
    #     else:
    #         return self.get_oriented_edge(first_interval, 
    #                 separatrix.traversed_intervals[1], 'center',
    #                 crossing = separatrix.intersections[0])

#xs    def _repr_(self):
        

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
            
        
                                  
                                
                                  


