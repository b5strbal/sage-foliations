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

from sage.structure.sage_object import SageObject
from train_track import TrainTrackPath
from base import LEFT

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
        self._tt = foliation.train_track()
        assert(not self._foliation.is_bottom_side_moebius() or 
                interval.side == 0) # otherwise it is not a real separatrix
        self._flip_count = 0
        interval = interval.add_to_position(end, foliation)
        p = interval.endpoint(LEFT, foliation)
        self._end_side = interval.endpoint_side(LEFT, foliation)

        self._intersections = [p]
        self._tt_path = TrainTrackPath()
        self._center_lengthen(p, interval)


        other_int = interval.prev(foliation)
        new_int = self._tt_path[-1].end()
        self._other_first_edge = self._tt.get_oriented_edge(other_int,
                                                            new_int,
                                                            'center',p,
                                                            foliation)

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
        # if self._foliation.is_bottom_side_moebius() and self._next_intersection != None:
        #     self._intersections.append(self._next_intersection)
        #     self._next_intersection = None
        #     return
        
        # moving to pair (and crossing center line if not moebius)
        last_int = self._tt_path[-1].end()
        if last_int.is_flipped(self._foliation):
            self._flip_count += 1
        self._end_side, p, new_int = self._foliation.apply_iet(self._intersections[-1],
                                                               self._end_side,
                                                               last_int)
        # print self._end_side
        self._intersections.append(p)
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                        'pair'))

        self._center_lengthen(p)

    def _center_lengthen(self, p, last_int = None):
        if last_int == None:
            last_int = self._tt_path[-1].end()

        self._end_side = (self._end_side + 1) % 2
        # print self._end_side, p
        # print self._intersections
        # print self._tt_path
        new_int = self._foliation.in_which_interval(p, self._end_side)
        self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
                                                        'center', p,
                                                        self._foliation))



        # p, new_int = self._foliation.point_int_on_other_side(p,
        #                                             last_int.side)
        # if self._foliation.is_bottom_side_moebius():
        #     self._next_intersection = p
        # else:
        #     self._intersections.append(p)
        # # print last_int, new_int
        # # print self._tt_path
        # # print self._intersections, self._next_intersection
        # self._tt_path.append(self._tt.get_oriented_edge(last_int, new_int,
        #                                             'center', p))
                
    def tt_path(self, end, endpoint = None, to_keep = 'before'):
        first_edge = self._tt_path[0] if end == 0 else \
                self._other_first_edge
        path = TrainTrackPath([first_edge] + self._tt_path[1:])
        if endpoint == None:
            return path
        # return path
        # print self._tt_path
        # print self._intersections
        # print endpoint

        cutting_index = 0
        # print endpoint
        while self._intersections[cutting_index] != endpoint:
            # print cutting_index, self._intersections[cutting_index]
            cutting_index += 1

        if to_keep == 'before':
            return TrainTrackPath(path[:2*cutting_index+1])
        if to_keep == 'after':
            return TrainTrackPath(path[2*cutting_index:])

    def _repr_(self):
        s = "Intersections: "
        s += repr(self._intersections)
        # s += "\nTrain Track Path: " + repr(self._tt_path)
        s += "\nTraversed intervals; "
        for oriented_edge in self._tt_path:
            s += repr(oriented_edge.start()) + ", "
        s += repr(self._tt_path[-1].end())
        s += "\nFlips: " + repr(self._flip_count)
        return s

    def intersections(self):
        return self._intersections

    # def intersections_without_endpoint(self):
    #     c = -1 if self._foliation.is_bottom_side_moebius() else -2
    #     return self._intersections[:c]

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
        # if the bottom side is moebius, an alternative way has to be considered
        # when trying to orient the surface by a double cover
        assert(not self._foliation.is_bottom_side_moebius())
        return (self.end_side == self.start_side)\
            == self.is_flipped()
        
    def raw_end_side(self):
        return self._tt_path[-1].start().side

    def end_side(self):
        return (self._end_side + 1) % 2

    def raw_start_side(self):
        return self._tt_path[0].start().side

    # @property
    # def start_point(self):
    #     return self._intersections[0]

    @property
    def endpoint(self):
        return self._intersections[-1]

    # def raw_endpoint(self):
    #     return self._foliation.adj_to_raw(self.endpoint,
    #                                       self.end_side())

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

