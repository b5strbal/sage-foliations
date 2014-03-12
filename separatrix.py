from sage.structure.sage_object import SageObject
from foliation import mod_one

class Separatrix(SageObject):
    def __init__(self, train_track, interval, end = 0, stop_function)
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
        Moreover, it is only 2 singularities, not 12, becuase vertices
        are identified. It is just that when the Foliation is drawn,
        then there is 12 problematic points, corresponding to the 12
        separatrices.

        The separatrices are in one-to-one correspondence with the 
        starting points of the intervals in the interval exchange.
        The side and position of a separatrix is defined to be the 
        side and position of the corresponding interval.

        """
        self._foliation = train_track.foliation
        self._tt = train_track
        assert(not self._foliation.is_bottom_side_moebius() or 
                interval.side == 0) # otherwise it is not a real separatrix
        self._flip_count = 0
        p = interval.endpoint(end)
        self._intersections = [p]
        self._tt_path = TrainTrack.Path()
        self._center_lengthen(p, interval)

        # this is -1 if we are at the last element of self._intersections
        # and -2 if it is the second to last element. Only used when the
        # bottom side is a Moebius band.
        self._current_index = -2

        while not stop_function():
            self.lengthen()


        
    def lengthen():
        if self._foliation.is_bottom_side_moebius() and current_index == -2:
            self._current_index = -1
            return
        
        self._current_index = -2
        last_int = self._tt_path[-1].endpoint(1)
        if last_int.is_flipped():
            self._flip_count += 1
        p, new_int = self._foliation.apply_iet(self._intersections[-1],
                                                    last_int)
        self._intersections.append(p)
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                        'pair'))
        self._center_lengthen(p)

    def _center_lengthen(self, p, last_int = self._tt_path[-1].endpoint(1)):
        p, new_int = self._foliation.point_int_on_other_side(p,
                                                    last_int.side)
        self._intersections.append(p)
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                    'center', p))
                

    @classmethod
    def from_bounding_arc(cls, train_track, interval, end, bounding_arc):
        def stop_function(self):
            return bounding_arc.contains(self._intersections[self._current_index])
        return Separatrix(train_track, interval, end, stop_function)

    @classmethod
    def from_flips_to_stop(cls, train_track, interval, end,
                           number_of_flips_to_stop):
        def stop_function(self):
            return self._flip_count == number_of_flips_to_stop
        return Separatrix(train_track, interval, end, stop_function)
        
    @classmethod
    def get_first_orientation_reversing(cls, train_track, interval, end):
        def stop_function(self):
            return self.is_orientation_reversing()
        return Separatrix(train_track, interval, end, stop_function)
    

    def _repr_(self):
        s = "Intersections: " + repr(self._intersections)
        s += "\nTraversed ints: " + repr(self._traversed_intervals)
        s += "\nFlips: " + repr(self._flip_count)
        return s

    @property
    def intersections(self):
        return self._intersections

    @property
    def foliation(self):
        return self._foliation

    def is_orientation_reversing(self):
        return (self.end_side == self.start_side)\
            == self.is_flipped()
        
    @property
    def end_side(self):
        if self.foliation.is_bottom_side_moebius() and \
           self._current_index == -1:
            return 1
        return self._tt_path[-1].start().side
    
    @property
    def start_side(self):
        return self._tt_path[0].start().side

    @property
    def endpoint(self):
        return self._intersections[self._current_index]

    # def first_interval(self):
    #     return self._traversed_intervals[0]

    def first_interval_end(self):
        return self._first_interval_end

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

    @classmethod
    def sorted_separatrices(cls, separatrices, starting_point = 0):
        return [sorted([s for s in separatrices if 
                        s.end_side == side],
                       key = lambda s: mod_one(s.endpoint() -
                                               starting_point))
                for side in {0, 1}]
        
        


    def shorten(self):
        self._traversed_intervals.pop()
        if self._traversed_intervals[-1].is_flipped():
            self._flip_count -= 1
        self._traversed_intervals.pop()
        self._intersections.pop()

    def shift_to_end(self, new_end):
        if new_end == self._first_interval_end:
            return
        if new_end == 1: #shift to left
            shift = -1
        if new_end == 0: #shift to right
            shift = 1
        self._first_interval_end = new_end
        self._traversed_intervals[0] =\
            self._traversed_intervals[0].add_to_position(shift)
        
    def _latex_(self):
        from foliation_latex import FoliationLatex
        return FoliationLatex(self._foliation).tikz_picture(
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
