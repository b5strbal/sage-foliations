from sage.structure.sage_object import SageObject
from mymath import mod_one
from train_track import TrainTrack
from constants import *

class Separatrix(SageObject):
    def __init__(self, foliation, interval, end = 0, bounding_arc = None,
                 number_of_flips_to_stop = None,
                 stop_at_first_orientation_reverse = False):
        r"""
        Calculates the segment of a separatrix before first intersecting
        a set of intervals.

        The separatrices considered here emanate from a singularity towards
        our distinguished simple closed curve, the interval $[0,1)$. For 
        example, a Foliation with Involution('a a b b c c') has six
        singularities "above the simple closed curve" and none below it.
        This terminology might be a little misleading, because 
        a Foliation with Involution('1 2 3 4 5 6', '2 1 4 3 6 5') has
        six singularities above and below as well, but here our curve is
        not separating, so the top and bottom part is actually connected.
        Moreover, it is only 2 singularities, not 12, because vertices
        are identified. It is just that when the Foliation is drawn,
        then there is 12 problematic points, corresponding to the 12
        separatrices.

        The separatrices are in one-to-one correspondence with the 
        starting points of the intervals in the interval exchange.
        The side and position of a separatrix is defined to be the 
        side and position of the corresponding interval.

        """
        self._foliation = foliation
        self._tt = foliation.train_track
        assert(not self._foliation.is_bottom_side_moebius() or 
                interval.side == 0) # otherwise it is not a real separatrix
        self._flip_count = 0
        interval = interval.add_to_position(end)
        p = interval.endpoint(0)
        self._intersections = [p]
        self._tt_path = TrainTrack.Path()
        self._next_intersection = None
        self._center_lengthen(p, interval)

        other_int = interval.prev()
        self._other_first_edge = self._tt.get_oriented_edge(other_int,
                                            self._tt_path[0].end(),
                                                             'center',p)
        
        
        # this is -1 if we are at the last element of self._intersections
        # and -2 if it is the second to last element. Only used when the
        # bottom side is a Moebius band.
        # self._current_index = -2

        if bounding_arc == None:
            if number_of_flips_to_stop != None:
                def terminate():
                    return self._flip_count == number_of_flips_to_stop
            elif stop_at_first_orientation_reverse:
                def terminate():
                    return self.is_orientation_reversing()
            else:
                assert(False)
        else:
            def terminate():
                return bounding_arc.contains(
                    self._intersections[-1])


        while not terminate():
            self.lengthen()



        
    def lengthen(self):
        # crossing the centerline
        if self._foliation.is_bottom_side_moebius() and self._next_intersection != None:
            self._intersections.append(self._next_intersection)
            self._next_intersection = None
            return
        
        # moving to pair (and crossing center line if not moebius)
        last_int = self._tt_path[-1].end()
        if last_int.is_flipped():
            self._flip_count += 1
        p, new_int = self._foliation.apply_iet(self._intersections[-1],
                                                    last_int)
        self._intersections.append(p)
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                        'pair'))
        self._center_lengthen(p)

    def _center_lengthen(self, p, last_int = None):
        if last_int == None:
            last_int = self._tt_path[-1].end()
        p, new_int = self._foliation.point_int_on_other_side(p,
                                                    last_int.side)
        if self._foliation.is_bottom_side_moebius():
            self._next_intersection = p
        else:
            self._intersections.append(p)
        # print last_int, new_int
        # print self._tt_path
        # print self._intersections, self._next_intersection
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                    'center', p))
                
    def tt_path(self, end, endpoint = None, to_keep = 'before'):
        first_edge = self._tt_path[0] if end == 0 else \
                self._other_first_edge
        path = TrainTrack.Path([first_edge] + self._tt_path[1:])
        if endpoint == None:
            return path
        # return path

        cutting_index = 0
        while self._intersections[cutting_index] != endpoint:
            cutting_index += 1

        if cutting_index % 2 == 1:
            cutting_index -= 1
        if to_keep == 'before':
            return TrainTrack.Path(path[:cutting_index+1])
        if to_keep == 'after':
            return TrainTrack.Path(path[cutting_index:])

    def _repr_(self):
        s = "Intersections: "
        s += repr(self.intersections())
        # s += "\nTrain Track Path: " + repr(self._tt_path)
        s += "\nTraversed intervals; "
        for oriented_edge in self._tt_path:
            s += repr(oriented_edge.start()) + ", "
        s += repr(self._tt_path[-1].end())
        s += "\nFlips: " + repr(self._flip_count)
        return s

    def intersections(self):
        return self._intersections

    def intersections_without_endpoint(self):
        c = -1 if self._foliation.is_bottom_side_moebius() else -2
        return self._intersections[:c]

    def get_intersection(self, n):
        return self._intersections[n]

    def get_tt_edge(self, n):
        return self._tt_path[n]

    def num_intersections(self):
        return len(self._intersections)

    @property
    def foliation(self):
        return self._foliation

    def is_orientation_reversing(self):
        return (self.end_side == self.start_side)\
            == self.is_flipped()
        
    @property
    def end_side(self):
        if self.foliation.is_bottom_side_moebius() and \
           self.num_intersections() % 2 == 0:
            return BOTTOM
        return self._tt_path[-1].start().side
    
    @property
    def start_side(self):
        return self._tt_path[0].start().side

    # @property
    # def start_point(self):
    #     return self._intersections[0]

    @property
    def endpoint(self):
        return self._intersections[-1]

    def first_edge(self, end):
        return self._tt_path[0] if end == 0 else self._other_first_edge

    def first_interval(self, end):
        return self.first_edge(end).start()

    # def first_interval_end(self):
    #     return self._end

    # def final_end(self):
    #     if not self.is_flipped():
    #         return self._end
    #     return (self._end + 1) % 2

    def is_flipped(self):
        return self._flip_count % 2 == 1

    @classmethod
    def get_all(cls, foliation, bounding_arc = None,
                number_of_flips_to_stop = None,
                stop_at_first_orientation_reverse = False):
        return [Separatrix(foliation, interval, 0, bounding_arc,
                           number_of_flips_to_stop,
                           stop_at_first_orientation_reverse) for
                interval in foliation.intervals()]


        
        

        
    def _latex_(self):
        return self._foliation.latex_options().tikz_picture(
            separatrices = [self])




    # def closing_intervals(self):
    #     repeats = 1
    #     if foliation.is_bottom_side_moebius() and self.end_side == 0:
    #         repeats = 2
    #     for i in range(repeats):
    #         self.lengthen()

    #     closing_intervals = self._traversed_intervals[-2*repeats : -1]
    #     del self._traversed_intervals[-2*repeats:]
    #     del self._intersections[-repeats:]
    #     return closing_intervals

    # def lengthened(self):
    #     from copy import copy
    #     new_separatrix = copy(self)
    #     new_separatrix._intersections = self._intersections[:]
    #     new_separatrix._traversed_intervals = self._traversed_intervals[:]
    #     new_separatrix.lengthen()
    #     return new_separatrix



    # def shorten(self):
    #     self._traversed_intervals.pop()
    #     if self._traversed_intervals[-1].is_flipped():
    #         self._flip_count -= 1
    #     self._traversed_intervals.pop()
    #     self._intersections.pop()

    # def shift_to_end(self, new_end):
    #     if new_end == self._first_interval_end:
    #         return
    #     if new_end == 1: #shift to left
    #         shift = -1
    #     if new_end == 0: #shift to right
    #         shift = 1
    #     self._first_interval_end = new_end
    #     self._traversed_intervals[0] =\
    #         self._traversed_intervals[0].add_to_position(shift)
