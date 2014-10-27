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



from separatrix import Separatrix
from arc import Arc
from collections import namedtuple
from sage.structure.sage_object import SageObject
from interval import Interval
from base import RestrictionError, SaddleConnectionError, LEFT, RIGHT


Coding = namedtuple("Coding", "side, index, end, num_flips1, num_flips2")


class TransverseCurve(SageObject):
    def __init__(self, foliation, coding):
        if coding.index >= foliation.num_intervals(coding.side):
            raise RestrictionError("Invalid interval index.")
        interval = Interval(coding.side, coding.index)

        # Checking if the curve is transverse
        if ((coding.num_flips1 % 2 == 1) != \
           (coding.num_flips2 % 2 == 1)) != \
                    interval.is_flipped(foliation):
            raise RestrictionError("Curve is not transverse to the foliation")


        other_int = interval.pair(foliation)
        other_end = coding.end
        if interval.is_flipped(foliation):
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
        # print self._sep
        arcs = [Arc(self._sep[i].endpoint, self._sep[(i+1)%2].endpoint,
                    *openness) for i in range(2)]
            

        self._arc = Arc(sep1.endpoint, sep2.endpoint, *openness)

        intersections = sep1.intersections()[:-1] + sep2.intersections()[:-1]

        # if one of the separatrices is longer, both arcs are considered,
        # so if the current one is wrong, we try the other one
        if coding.num_flips1 > 0 or coding.num_flips2 > 0:
            x = intersections[0]
            if self._arc.contains(x):
                self._arc = Arc(sep2.endpoint, sep1.endpoint, *openness)
                self._sep = [sep2, sep1]

        # if self._arc.length() > 0.5 and sep1.end_side == sep2.end_side:
        #     raise RestrictionError("There is either no such curve without "
        #                            "self-intersection or "
        #                            "the curve is one-sided and has length"
        #                            " greater than half, so we can't compute"
        #                            " the train track transition map. ")
        # print self._arc
        # print sep1, sep2
        for x in intersections:
            # print x
            if self._arc.contains(x):
                raise RestrictionError("The curve is self-intersecting")

        self._direction = RIGHT if openness[0] == 'closed' else LEFT
        self._coding = coding
        # print coding, self._direction, self._arc
        # print '\n'

    def _repr_(self):
        s = ''
        s += 'Arc: ' + repr(self._arc) #+ '(' + repr(self._arc.length()) + ')'
        return s

    def _latex_(self):
        return self._foliation.latex_options().tikz_picture(
            transverse_curves = [self], separatrices = self._get_separatrices())

    def _get_separatrices(self):
        return Separatrix.get_all(self._foliation, self._arc)

    @property
    def _foliation(self):
        return self._sep[0].foliation

    def separatrix(self, n):
        return self._sep[n]

    def direction(self):
        return self._direction

    def arc(self):
        return self._arc

    def coding(self):
        return self._coding
                           
    def is_one_sided(self):
        is_same_side = self._sep[0].end_side() == self._sep[1].end_side()
        return self.is_twisted() != is_same_side

    def is_twisted(self):
        return self._foliation.is_bottom_side_moebius() and \
           self._arc[0] > self._arc[1]

        
        

    def new_foliation(self):
        from transition_map import new_foliation
        fol = self._foliation
        adj_starting_side = self._sep[0].end_side() if self._direction == RIGHT \
                       else self._sep[1].end_side()
        # self.save('tc.save')
        return new_foliation(self._get_separatrices(), self._sep[0].endpoint,
                             adj_starting_side,
                             is_one_sided = self.is_one_sided(),
                             hdir = self._direction,
                             adj_ending_point = self._sep[1].endpoint,
                             transformation_coding = RestrictionCoding(fol, self._coding))


RestrictionCoding = namedtuple("RestrictionCoding", "foliation, coding")


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



def other_tts(fol, depth):
    if depth <= 0:
        if fol.train_track()._degrees == [3,3,3,3,3,3]:
            print fol
        return
    codings = get_codings(fol)
    for c in codings:
        tc = TransverseCurve(fol, c)
        other_tts(tc.new_foliation()[0], depth-1)
