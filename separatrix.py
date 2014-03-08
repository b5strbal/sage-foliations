from sage.structure.sage_object import SageObject
from foliation import mod_one

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
        This terminology might be a little misleading, beacuse 
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

        INPUT:

        - ``side`` - 0 or 1, the side of the separatrix, 
            the 0 for top, 1 for bottom side

        - ``pos`` - non-negative integer, the position of the separatrix

        - ``intervals`` - list of Interval objects that stand for "walls"
            here, i.e. the separatrix is followed until it hits of
            there (potentially overlapping) walls. The default value,
            None carries a different meaning: in this case the process
            stops when the first flipped interval is encountered, i.e.
            when the orientation is first reversed on the separatrix

        OUTPUT:

        - Separatrix -- it is a namedtuple with two attributes:
            "intersections" and "is_flipped". "intersections" is a list
            of interval positions of intersections and points of 
            intersections (see examples). "is_flipped" is a boolean,
            which says if the separatrix segment is flipped, i.e. 
            reverses orientation (if the Intervals argument is None,
            this is always the case by definition)

        EXAMPLES::

            sage: from sage.dynamics.foliations.foliation import Interval
            sage: f = Foliation.nonorientable_arnoux_yoccoz(4)
            sage: i = Interval(f._divpoints[0][0], f._divpoints[0][2])
            sage: f._first_intersection(0, 0, [i])
            Separatrix(intersections=[0], traversed_intervals=[(0, 0)], 
                is_flipped=False)
            sage: f._first_intersection(0, 1, [i])
            Separatrix(intersections=[0.271844506346], 
                traversed_intervals=[(0, 1)], is_flipped=False)
            sage: f._first_intersection(0, 2, [i])
            Separatrix(intersections=[0.543689012692, 0.0436890126921], 
                traversed_intervals=[(0, 2), (1, 1), (1, 0)], is_flipped=False)

            sage: f = Foliation.RP2_arnoux_yoccoz()
            sage: f._first_intersection(0, 0)
            Separatrix(intersections=[0, 0.5, 0.567442248868], 
                traversed_intervals=[(0, 0), (1, 0), (1, 1), (0, 2), (0, 3)], 
                is_flipped=True)

        """
        assert(not foliation.is_bottom_side_moebius() or 
                interval.side == 0) # otherwise it is not a real separatrix
        self._flip_count = 0
        self._intersections = [interval.endpoint(end)] 
        self._traversed_intervals = [interval]
        self._foliation = foliation
        self._first_interval_end = end

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
                return bounding_arc.contains(self._intersections[-1])

        while not terminate():
            self._lengthen()

    def _repr_(self):
        s = "Intersections: " + repr(self._intersections)
        s += "\nTraversed ints: " + repr(self._traversed_intervals)
        s += "\nFlips: " + repr(self._flip_count)
        return s

    @property
    def intersections(self):
        return self._intersections

    @property
    def traversed_intervals(self):
        return self._traversed_intervals

    @property
    def foliation(self):
        return self._foliation

    def is_orientation_reversing(self):
        return (self.end_side == self._traversed_intervals[0].side)\
            == self.is_flipped()
        
    @property
    def end_side(self):
        return self.traversed_intervals[-1].side
    
    def endpoint(self):
        return self._intersections[-1]

    def first_interval(self):
        return self._traversed_intervals[0]

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
        
        
        

    def _interval_on_other_side(self):
        return self._foliation.in_which_interval(\
                self._intersections[-1], (self.end_side + 1) % 2)
                
    def _lengthen(self):
        self._traversed_intervals.append(self._interval_on_other_side())
        if self._traversed_intervals[-1].is_flipped():
            self._flip_count += 1
        point, interval = self._foliation.apply_iet(\
                self._intersections[-1], 
                self._traversed_intervals[-1])
        self._traversed_intervals.append(interval)
        self._intersections.append(point)

    def closing_intervals(self):
        repeats = 1
        if foliation.is_bottom_side_moebius() and self.end_side == 1:
            repeats = 2
        for i in range(repeats):
            self._lengthen()

        closing_intervals = self._traversed_intervals[-2*repeats : -1]
        del self._traversed_intervals[-2*repeats:]
        del self._intersections[-repeats:]
        return closing_intervals

    def lengthened(self):
        from copy import copy
        new_separatrix = copy(self)
        new_separatrix._intersections = self._intersections[:]
        new_separatrix._traversed_intervals = self._traversed_intervals[:]
        new_separatrix._lengthen()
        return new_separatrix

    def _latex_(self):
        from foliation_latex import FoliationLatex
        return FoliationLatex(self._foliation).tikz_picture(
            separatrices = [self])
