
from sage.graphs.digraph import DiGraph
from foliation import mod_one

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

    class OrientedEdge(namedtuple('OrientedEdge', 'edge, direction')):
        def __new__(cls, edge, direction):
            return super(OrientedEdge, cls).__new__(cls, edge, direction)

        def reversed(self):
            return OrientedEdge(self.edge, -self.direction)

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

    def get_first_edge(self, separatrix, left_or_right):
        first_interval = separatrix.traversed_intervals[0]
        if left_or_right == 'left':
            first_interval = first_interval.prev()
        if self._foliation.is_bottom_side_moebius():
            return self.get_oriented_edge(first_interval, 
                    separatrix.traversed_intervals[3], 'center')
        else:
            return self.get_oriented_edge(first_interval, 
                    separatrix.traversed_intervals[1], 'center',
                    crossing = separatrix.intersections[0])

    class Path:
        def __init__(self, oriented_edges):
            self._path = oriented_edges

        @classmethod
        def from_separatrix(cls, separatrix):
            closing_ints = separatrix.closing_intervals()
            trav_ints = separatrix.traversed_intervals
            trav_ints.extend(closing_ints)
            path = []
            if separatrix._foliation.is_bottom_side_moebius():
                for i in range(1, (len(separatrix.intersections) + 1) // 2):
                    path.append(self.get_oriented_edge(trav_ints[4 * i - 1],
                        trav_ints[4 * i], 'pair'))
                    path.append(self.get_oriented_edge(trav_ints[4 * i],
                        trav_ints[4 * i + 3], 'center'))
            else:
                for i in range(1, len(separatrix.intersections)):
                    path.append(self.get_oriented_edge(trav_ints[2 * i - 1], 
                        trav_ints[2 * i], 'pair'))
                    path.append(self.get_oriented_edge(trav_ints[2 * i],
                        trav_ints[2 * i + 1], 'center', 
                        crossing = separatrix.intersections[i]))
            del trav_ints[-len(closing_ints):]
            return cls(path)

        def reversed(self):
            return TrainTrack.Path([oriented_edge.reversed() for 
                oriented_edge in reversed(self._path)])

        def __add__(self, other):
            assert(self._path[-1] == other._path[0])
            return TrainTrack.Path(self._path[:-1] + other._path)






