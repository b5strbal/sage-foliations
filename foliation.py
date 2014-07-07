r"""


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

from sage.dynamics.interval_exchanges.constructors import GeneralizedPermutation
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import vector
from base import *
import random

from interval import Interval
from sage.rings.rational import Rational
from collections import namedtuple

SymmetryCoding = namedtuple("SymmetryCoding", "interval, hdir")


class Foliation(SageObject):
    r"""
    A singular measured foliation on a surface.

    Given a singular measured foliation we consider a simple closed
    curve `\gamma` on the surface which is transverse to the foliation. A small
    neighborhood of `\gamma` is a Moebius strip or an annulus
    depending on whether `\gamma` is one-sided or two-sided. The
    foliation defines an interval exchange on the boundary `B` of this
    neighborhood defined as follows. Consider the first intersections
    of the separatrices of the foliation with `B`. They divide the
    circle componenents of `B` to intervals, which are identified by
    the foliation is pairs. 

    NOTE: We will assume that there is at least one division point
    on each component of `B`. Therefore foliations on the closed torus
    or the Klein bottle -- that don't have singularities -- cannot be
    represented. In general, foliations on surfaces with negative Euler
    characteristic (including the punctured tori and Klein bottles)
    are considered.

    NOTE: Another convention is that the surface is assumed to have as
    few punctures are allowed by the foliation, i.e. punctures are
    exactly at 1-pronged and 2-pronged singularities. (1-pronged
    singularities are only allowed at punctures, and a 2-pronged
    singularity is only a singularity if it is at a puncture.) At this
    point there is no way to represent a foliation that has a
    3-pronged singularity at a puncture, for example. 

    We normalize the measure such that components of `B` have measure
    1. So the measure of `\gamma` is 1 or 1/2 if it is two-sided or
    one-sided, respectively. This way we can parametrize the division
    points with numbers between 0 and 1. We choose 0 to be one of the
    division points. There are as many choices as division
    points. There are also two ways of choosing the direction of the
    identification of `[0,1]` with components of `B`. 

    If `B` has two components, we will call one the top side, the
    other the bottom side. When there is one component, we will call
    that the top side, and say that there is a Moebius band on the
    bottom side. Label the intervals on each side so that pairs have
    the same label. Some pairs of intervals are identified by a
    translation, some by a flips. We will define the foliation by the
    list of top and bottom labels, the set of flipped labels, the
    lengths of the intervals, and the twist in case `B` has two
    components.

    The conventions for defining and printing interval exchanges are
    the same as in iet.GeneralizedPermutation.

    INPUT:

    - ``top_labels``, ``bottom_letters`` -- either strings of the
      labels separated by spaces or lists of the label strings. The
      two list correspond to the two components of `B`. To indicate
      that there is one component, i.e. that `\gamma` is one-sided,
      set ``bottom_letters`` to ``moebius``. 


    - ``lengths`` -- (default: None) as in
    iet.IntervalExchangeTransformation, this is either a list, tuple
    or dictionary of the length parameters.  If it is a list or
    tuple, then the lengths are assigned to intervals in the order
    of their appearance in the ``top_letters`` and
    ``bottom_letters``. The length of the list or tuple should equal
    the number of labels. If it's a dict, it should contain entries
    in the form label:length. Finally, if omitted, lengths are
    generated randomly.

    - ``flips`` -- (default: []) the list of flipped labels. It the
      labels consist of single characters, ``flips`` can be a string
      containing those characters.
    
    - ``twist`` -- (default:None) a real number, the twist
    parameter. When `\gamma` is one-sided, this argument is ignored,
    because `B` has only one component, so the twist doesn't have a
    meaning. Otherwise it should be some non-zero number (otherwise
    there is a saddle connection). Omitting it generates a random twist.

    EXAMPLES:
    
    Here are two different but equivalent definitions::

        sage: f = Foliation('a a b b', 'c c', [1, 2, 3], twist=1/2, flips = 'b'); f
         a a -b -b
         c c
        Lengths: {'1': 1/6, '3': 1/2, '2': 1/3}
        Twist: 1/12

        sage: g = Foliation(['a','a','b','b'], ['c','c'], {'a':1,'b':2,'c':3}, twist=1/2, flips = 'b'); g
         a a -b -b
         c c
        Lengths: {'1': 1/6, '3': 1/2, '2': 1/3}
        Twist: 1/12

        sage: f == g
        True

    Omitting lengths and/or twist is fine, they are generated randomly::

        sage: f = Foliation('a a b b', 'c c')

    If `\gamma` is two-sided, the lengths of the two components of `B`
    are equal, so in that sense certain length specifications are
    illegal. Instad of raising an error, however, if the sum of the
    lengths on top and bottom are not equal, the lengths are
    automatically adjusted to make it work::

        sage: Foliation('a a b b','c c', [0.2, 0.3, 0.5000001], twist
        = 0.15)
        a a b b
        c c
        Lengths: {'a': 0.200000059999988, 'c': 0.500000000000000, 'b': 0.299999940000012}
        Twist: 0.149999970000006
      
    For easier comparison of two ``Foliation`` objects, the twist is
    normalized to a positive number as small as possible by cyclically
    permuting the bottom labels::

       sage: f=Foliation('a b c','c b a', [0.1, 0.2, 0.7], twist =
       0.25); f
       a b c
       a c b
       Lengths: {'a': 0.100000000000000, 'c': 0.700000000000000, 'b': 0.200000000000000}
       Twist: 0.150000000000000
    
    Saddle connections are not checked in the constructor, for
    instance the following doesn't give an error::

       sage: f=Foliation('a b c d','d c b a', [0.25, 0.25, 0.25, 0.25], twist = 0.25)

    However, certain methods will raise an error when they realize the
    saddle connections. Certain permutations result in saddle
    connections for arbitrary length and twist parameters. These
    permutations are not straightforward to detect, so no error is
    given. For instance 'a a b b', 'c c' is such a permutation.

    """ 

    def __init__(self, top_labels, bottom_letters, lengths = None, flips = [], 
            twist = None):
        """
        TESTS::
            sage: f = Foliation('a b c','a c b', [1, 2, 3], twist=1); f
            a b c
            a c b
            Lengths: {'a': 1/6, 'c': 1/2, 'b': 1/3}
            Twist: 1/6
            sage: f._divpoints
            [[0, 1/6, 1/2], [1/6, 1/3, 5/6]]
            sage: f._lengths['a']
            1/6
            sage: f._lengths['b']
            1/3
            sage: f._lengths['c']
            1/2
            sage: f._twist
            1/6
            sage: f._gen_perm
            a b c
            a c b
            sage: f._is_curve_one_sided
            False
        
            sage: f = Foliation('a a b b c c','moebius', [1, 2, 3]); f
            a a b b c c
            Moebius band
            Lengths: {'a': 1/12, 'c': 1/4, 'b': 1/6}
            sage: f._divpoints
            [[0, 1/12, 1/6, 1/3, 1/2, 3/4], [0, 1/2]]
            sage: f._lengths['a']
            1/12
            sage: f._lengths['b']
            1/6
            sage: f._lengths['c']
            1/4
            sage: f._twist
            0
            sage: f._gen_perm
            a a b b c c
            JOKER JOKER
            sage: f._is_curve_one_sided
            True

        """
        from bisect import bisect_left

        self._is_bottom_side_moebius = False
        if bottom_letters == 'moebius': # bottom side is Moebius
            bottom_letters = 'JOKER JOKER'
            self._is_bottom_side_moebius = True

        # initializing the permutation, self._gen_perm
        self._init_gen_perm(top_labels, bottom_letters, flips)

        # initializing the lengths, self._lengths
        if lengths is None:
            while True:
                lengths = [random.uniform(EPSILON, 1.0)
                           for i in range(len(self._labels))]
                try:
                    total = self._init_lengths(lengths)
                    break
                except ValueError:
                    pass
        else:
            total = self._init_lengths(lengths)

        # setting/normalizing twist
        if twist is None:
            twist = random.uniform(EPSILON, total)
        if self.is_bottom_side_moebius():
            twist = 0
        twist /= total

        # setting top divpoints, self._divpoints, the bottom is not
        # yet twisted
        self._divpoints = [[0], [0]]
        for interval in self._all_intervals[TOP] + self._all_intervals[BOTTOM]:
            self._divpoints[interval.side].append(self._divpoints[
                interval.side][-1] + self._lengths[interval.label(self)])
        for side in [TOP,BOTTOM]:
            self._divpoints[side].pop()

        preimage_of_zero = mod_one(-twist)
        containing_int = bisect_left(self._divpoints[BOTTOM], 
                                     preimage_of_zero) % len(self._divpoints[BOTTOM])

        # reseting the permutation after normalizing the twist
        self._init_gen_perm(self._gen_perm_list[TOP],
                            self._gen_perm_list[BOTTOM][containing_int:]+
                            self._gen_perm_list[BOTTOM][:containing_int], 
                            flips = self._flips)
        self._twist = mod_one(self._divpoints[BOTTOM][containing_int] -
                preimage_of_zero)

        # reseting the bottom divpoints now with the twist
        self._divpoints[BOTTOM] = [self._twist]
        for interval in self._all_intervals[BOTTOM][:-1]:
            self._divpoints[1].append(self._divpoints[BOTTOM][-1] + 
                                      interval.length(self))


    def _init_gen_perm(self, top_labels, bottom_letters, flips):
        """Initialize the permutation-related variables.

        More precisely: _gen_perm, _gen_perm_list, _all_intervals,
        _numerical_label, _pair. Other initializing methods depend on
        these variables, therefore this one must be called first.

        INPUT:

        - ``top_labels``,``bottom_labels``,``flips`` -- as in the constructor

        """
        self._gen_perm = GeneralizedPermutation(\
                        top_labels, bottom_letters, flips = flips)

        # The permutation in list form. Access is faster.
        self._gen_perm_list = self._gen_perm.list()

        # Saving the set of labels.
        self._labels = set(self._gen_perm.alphabet())
        self._labels.discard('JOKER')
        
        # The list of flipped labels.
        self._flips = self._gen_perm.flips() if \
                      isinstance(self._gen_perm[0][0], tuple) else []
        # In the first case, self._gen_perm is a FlippedLabelledPermutationLI
        # that has a flips() method
        # In the second case, self._gen_perm is a LabelledPermutationIET 
        # that doesn't have a flips method, but then there are no
        # flips, so we assign the empty list.

        
        # initializing self._numerical_label and self._pair
        self._all_intervals = [[], []]
        label_to_interval = {}
        self._numerical_label = {}
        count = 0
        self._pair = {}
        for side in {0, 1}:
            for index in range(len(self._gen_perm_list[side])):
                interval = Interval(side, index)
                self._all_intervals[side].append(interval)
                label = interval.label(self)
                if label not in self._numerical_label:
                    self._numerical_label[label] = count
                    count += 1
                    label_to_interval[label] = interval
                else:
                    self._pair[label_to_interval[label]] = interval
                    self._pair[interval] = label_to_interval[label]
            

    def _init_lengths(self, lengths):
        """
        Initializes self._lengths.

        Should only be used in the constructor of Foliation.

        INPUT: 

        - ``lengths`` -- either a list, tuple or dict or lengths. The
          lengths have to be positive, but even then, after balancing
          the lengths out, some lengths might turn into negative, in
          which case a ValueError is raised.

        """
        if isinstance(lengths, (list, tuple)):
            if len(lengths) != len(self._labels):
                raise ValueError('Bad number of lengths')
            self._lengths = {label: lengths[self._numerical_label[label]] 
                    for label in self._labels}
        elif isinstance(lengths, dict):
            if set(self._labels) != set(lengths.keys()):
                raise ValueError('Invalid length specification')     
            self._lengths = dict(lengths)
        else:
            TypeError('The ``lengths`` argument should be a list, tuple'
                       'or dict')
            
        if self.is_bottom_side_moebius():
            self._lengths['JOKER'] = sum(self._lengths.values())

        totals = [sum(interval.length(self) for interval in 
                      self._all_intervals[side]) for side in [0,1]]

        #adjusting lengths in case they 
        #differ on top/bottom
        for interval in self._all_intervals[0]:
            if interval.pair(self).side == interval.side:
                self._lengths[interval.label(self)] += \
                        (totals[1] - totals[0])/2
                break

        if any(v <= 0 for v in self._lengths.values()):            
            raise ValueError('Corrected lengths are not positive')
        
        for label in self._lengths:
            self._lengths[label] /= totals[1]

        return totals[1]


    def _init_singularity_partition(self):
        """Initialize self._singularity_partition. """
        done = set()
        partition = []
        for interval in self.intervals():
            if interval in done:
                continue
            new_interval = interval
            partition.append([])
            end = 0 # we start from the left side of the interval
            while True:
                sing_index = new_interval if end == 0 else new_interval.next(self)
                if sing_index == interval and len(partition[-1]) > 0:
                    break
                done.add(sing_index)
                partition[-1].append(sing_index)
                end = (end + 1) % 2
                new_interval = new_interval.add_to_position((-1)**end, self)
                new_interval = new_interval.pair(self)
                if new_interval.is_flipped(self):
                    end = (end + 1) % 2
                
        self._singularity_partition = partition

    
    def permutation(self):
        """Return the permutation of the interval exchange.

        OUTPUT:

        a iet.GeneralizedPermutation object

        EXAMPLES::

        sage: f = Foliation('a a b b', 'c c', [1, 2, 3])
        sage: f.permutation()
        a a b b
        c c

        sage: f = Foliation('a a b b c c','moebius',flips='abc')
        sage: f.permutation()
        -a -a -b -b -c -c
         JOKER  JOKER
        
        """
        return self._gen_perm

    def intervals(self):
        """Return the list of intervals.

        This method is useful for iterating through the intervals of
        the foliation. The intervals are repserented as Interval
        objects with a range of convenient properties.

        OUTPUT:

        a list of ``Interval`` objects ordered as top intervals first,
        from left to right, then (if the bottom side is not a Moebius
        band) bottom intervals from left to
        right.

        EXAMPLES::

        sage: f = Foliation('a a b b', 'c c', [1, 2, 3])
        sage: f.intervals()
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1)]

        sage: f = Foliation('a a b b c c','moebius',flips='abc')
        sage: f.intervals()
        [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
        
        """
        if self.is_bottom_side_moebius():
            return self._all_intervals[0]
        return self._all_intervals[0] + self._all_intervals[1]


    def is_orientable(self, foliation_or_surface):

        """
        Test orientability of the foliation or surface.

        INPUT:

        - ``foliation_or_surface`` -- either 'foliation' or 'surface',
          depending on the orientability of the foliation or the
          surface is tested

        OUTPUT:

        True or False

        EXAMPLES:

        The foliation is orientable iff there are no flips::

            sage: f = Foliation('a a b b', 'c c', [1,2,3], flips ='a',
            twist = 1)
            sage: f.is_orientable('foliation')
            False

            sage: g = Foliation('a a b b', 'c c', [1,2,3], twist = 1)
            sage: g.is_orientable('foliation')
            True


        If the bottom side is a Moebius band, the surface is automatically
        non-orientable::

            sage: f = Foliation('a a b b','moebius',[1,2])
            sage: f.is_orientable('surface')
            False

        Also if there is flipped pair on different sides::

            sage: f = Foliation('a b c','c b a',[1,2,3], flips = 'a', twist=1)
            sage: f.is_orientable('surface')
            False

        Or if there is not flipped pair on the same side::

            sage: f = Foliation('a a b b', 'c c', [1,2,3],flips='ac',twist=1)
            sage: f.is_orientable('surface')
            False

        Otherwise it is orientable::

            sage: f = Foliation('a b c', 'b c a',[1,2,3],twist=1)
            sage: f.is_orientable('surface')
            True

        """
        if foliation_or_surface == 'foliation':
            return len(self._flips) == 0
        if foliation_or_surface == 'surface':
            # including the JOKER JOKER pair if the bottom side is Moebius 
            return all(not interval.is_orientation_reversing(self)
                       for x in self._all_intervals for interval in x)


    def num_intervals(self, side):
        """Return the number of intervals on the specified side.

        INPUT:

        - ``side`` -- (0=TOP or 1=BOTTOM) 

        OUTPUT:

        the number of intervals on the specified side.

        EXAMPLES::

        sage: f = Foliation('a a b b', 'c c', [1, 2, 3])
        sage: f.num_intervals(TOP)
        4
        sage: f.num_intervals(BOTTOM)
        2

        sage: f = Foliation('a a b b c c','moebius',flips='abc')
        sage: f.num_intervals(TOP)
        6
        sage: f.num_intervals(BOTTOM)
        0

        """
        if self.is_bottom_side_moebius() and side == BOTTOM:
            return 0
        return len(self._all_intervals[side])

    def train_track(self):
        """Return the canonical train track carrying the foliation.

        OUTPUT:

        The ``TrainTrack`` carrying the foliation.  If there is a
        saddle connection, the train track cannot be constructed, and
        SaddleConnectionError is raised.

        EXAMPLES::

        sage: f = Foliation('a a b b', 'c c', [1, 2, 3])
        sage: tt = f.train_track()

        """
        from train_track import TrainTrack
        if not hasattr(self, '_tt'):
            self._tt = TrainTrack(self)
        return self._tt


    def is_bottom_side_moebius(self):
        """
        Decides if the bottom side is Moebius band without
        punctures.

        This happens exactly when the second argument in
        the constructor is 'moebius', or in other words when the
        transverse curve `\gamma` is one-sided.

        OUTPUT:

        True or False

        EXAMPLES::

            sage: f = Foliation('a a b b','moebius',[1,2])
            sage: f.is_bottom_side_moebius()
            True

        Here the bottom component is a once punctured Moebius
        band, so it returns False::

            sage: f = Foliation('a a b b', 'c c',[1,2,3],twist=1)
            sage: f.is_bottom_side_moebius()
            False
        """
        return self._is_bottom_side_moebius


    def __eq__(self, other):
        r"""
        Tests exact equality.

        TESTS::

        sage: f = Foliation('a a b b c c','moebius',[1,2,3])
        sage: g = Foliation('1 1 2 2 3 3','moebius',[2.0000001,4,6])
        sage: h = Foliation('1 1 2 2 3 3','moebius',[2,4,6])
        sage: f == g
        False
        sage: f == h
        True
        
        """
        return self.equals(other,0)
                
    
    def equals(self, other, allowed_error = EPSILON):
        """Test approximate equality of two foliations.

        For the two foliations to be considered equal, the permutation
        of intervals has to be the same (but possibly with different
        labels). Also the vectors made from the length and twists
        parameters must differ by less than ``allowed_error``.

        INPUT:

        - ``other`` -- another ``Foliation`` object

        - ``allowed_error`` -- a positive number, the distance between
        the two length-twist vectors must be smaller than this in
        order for equality.
        
        EXAMPLES::

        sage: f = Foliation('a a b b c c','moebius',[1,2,3])
        sage: g = Foliation('1 1 2 2 3 3','moebius',[2.0000001,4,6])
        sage: f.equals(g,0.0001)
        True

        """
        return self.is_bottom_side_moebius() == other.is_bottom_side_moebius()\
                and self._gen_perm == other._gen_perm and \
                abs(vector(flatten(self._divpoints)) -
                vector(flatten(other._divpoints))) <= allowed_error
        
    def _repr_(self):
        r"""
        Return a representation of self.

        TESTS::

        sage: Foliation('a b c', 'a c b', [1, 2, 3], twist=1)
        a b c
        a c b
        Lengths: {'a': 1/6, 'c': 1/2, 'b': 1/3}
        Twist: 1/6

        sage: Foliation('a a b b c c','moebius', [2, 2, 1])
        a a b b c c
        Moebius band
        Lengths: {'a': 1/5, 'c': 1/10, 'b': 1/5}

        """
        if self.is_bottom_side_moebius():
            d = dict(self._lengths)
            del d['JOKER']
            return "{0}\nLengths: {1}".format(\
                repr(self._gen_perm).split('\n')[0] + '\nMoebius band', d)
        return "{0}\nLengths: {1}\nTwist: {2}".format(repr(self._gen_perm),
                self._lengths, self._twist)

    def _latex_(self):
        """Return a TikZ representation of the Foliation.

        OUTPUT:

        string, the TikZ representation of the Foliation.

        EXAMPLES::

        sage: f = Foliation('a a b b c c','moebius', [0.12312,0.531354,0.5134])
        sage: latex(f) #indirect doctest
        '...tikzpicture...'

        """
        return self.latex_options().tikz_picture()

    def latex_options(self):
        """Return the LaTeX options for the TikZ representation of
        self.

        OUTPUT:

        a ``FoliationLatex`` object, storing the set LaTeX options.

        EXAMPLES:

        Only non-default options are listed, so initially nothing is printed::
        
        sage: f = Foliation('a a b b c c','moebius', [0.12312,0.531354,0.5134])
        sage: f.latex_options()
        {}

        One can custom drawing options::

        sage: f = Foliation('a a b b c c','moebius', [0.12312,0.531354,0.5134])
        sage: f.set_latex_options(scale=10, interval_labelling=False)
        sage: f.latex_options()
        {'scale': 10, 'interval_labelling': False}

        """
        if not hasattr(self, '_latex_opts'):
            from foliation_latex import FoliationLatex
            self._latex_opts = FoliationLatex(self)
        return self._latex_opts

    def set_latex_options(self, **kwds):
        """Set LaTeX options.
        
        INPUT:

        - ``**kwds`` -- options to be set

        EXAMPLES::

        sage: f = Foliation('a a b b c c','moebius', [0.12312,0.531354,0.5134])
        sage: f.set_latex_options(scale=10, interval_labelling=False)
        sage: f.latex_options()
        {'scale': 10, 'interval_labelling': False}

        """
        opts = self.latex_options()
        opts.set_options(**kwds)

    def singularity_partition(self):
        """
        Return the singularity partition of the ``Foliation``.

        OUTPUT:

        - list of lists - each element of the list is a list corresponding to
          a singularity, and each such list contains the
          ``Interval``s whose left endpoint lies in a separatrix
          emanating from that singularity.

        EXAMPLES::

            sage: f = Foliation('a a b b c c','moebius')
            sage: sp = f.singularity_partition()
            sage: len(sp)
            1
            sage: set(sp[0]) == {(0,2), (0,3), (0, 4), (0, 5), (0, 0), (0, 1)}
            True

            sage: f = Foliation('a b c d', 'b a d c', flips = 'ac')
            sage: sp = f.singularity_partition()
            sage: len(sp)
            2
            sage: set(sp[0]) == {(1, 1), (0, 2), (1, 0), (0, 1)}
            True
            sage: set(sp[1]) == {(0, 0), (1, 3), (0, 3), (1, 2)}
            True

        """
        if not hasattr(self, '_singularity_partition'):
            self._init_singularity_partition()

        return list(self._singularity_partition)

    def singularity_type(self):
        """
        Return the singularity type of the ``Foliation``.

        OUTPUT:

        a tuple of the number of prongs of the singularities, in
        decreasing order.

        EXAMPLES::

            sage: f = Foliation('a a b b c c', flips = 'abc')
            sage: f.singularity_type()
            (3, 1, 1, 1)

            sage: f = Foliation('a a b b c c', 'd d', flips = 'abc')
            sage: f.singularity_type()
            (3, 2, 1, 1, 1)

        Finally, a once and twice punctured torus::

            sage: f = Foliation('a', 'a')
            sage: f.singularity_type()
            (2,)

            sage: f = Foliation('a b', 'a b')
            sage: f.singularity_type()
            (2, 2)

        """
        t = sorted([x for x in map(len, self.singularity_partition())],
                   reverse = True)
        return tuple(t)

    def euler_char(self):
        """Return the Euler characteristic of the closed surface.
        
        OUTPUT:

        an integer, the Euler characteristic

        EXAMPLES::

            sage: f = Foliation('a a b b c c','moebius',[0.12,0.45,0.63],flips = 'abc')
            sage: f.euler_char()
            -2

            sage: f = Foliation('a a b b c c', 'd d', flips = 'abc')
            sage: f.singularity_type()
            (3, 2, 1, 1, 1)

        A once punctured torus::

            sage: f = Foliation('a', 'a')
            sage: f.euler_char()
            0
        """
        ec = 0
        for p in self.singularity_type():
            ec += 2-p
        return ec/2

    def genus(self):
        """Return the genus of the surface.

        In the orientable case the genus of the number of tori in a
        connected sum representation, while in the non-orientable case
        it is the number of projective planes.

        OUTPUT:

        a non-negative integer, the genus of the surface

        EXAMPLES::

        """
        ec = self.euler_char()
        if self.is_orientable('surface'):
            return 1 - ec/2
        return 2 - ec

    def num_punctures(self):
        """Return the number of punctures.

        By convention, the surface is considered punctured at the 1-
        and 2-pronged singularities, but not at the others.

        OUTPUT:
        
        a non-negative integer, the number of punctures

        EXAMPLES::
        
        """
        return len([x for x in self.singularity_type() if x <= 2])

    def surface_type(self):
        """Return the surface type.

        Let `S_{g,n}` and `N_{g,n}` denote the orientable and
        nonorientable surface with genus `g` and `n` punctures,
        respectively. In the case `n=0` we will simply say `S_g` and `N_g`.

        OUTPUT:

        a string, the surface type

        EXAMPLES::

        """
        stem = 'S_' if self.is_orientable('surface') else 'N_'
        g = self.genus()
        n = self.num_punctures()
        if n == 0:
            return stem + str(g)
        return '{0}{{{g},{n}}}'.format(stem, g=g, n=n)


    
    def with_changed_lengths(self, length_vector = None):

        """Define a new Foliation with different length parameters.

        INPUT:

        - ``length_vector`` -- 


        """
        if self.is_bottom_side_moebius():
            return Foliation(self._gen_perm_list[0], 'moebius', length_vector,
                             flips = self._flips)
        if length_vector is None:
            lengths = None
            twist = None
        else:
            lengths = length_vector[:-1]
            twist = length_vector[-1]
        return Foliation(self._gen_perm_list[0], self._gen_perm_list[1],
                         lengths, flips = self._flips,
                         twist = twist)



    
    def raw_to_adj(self, raw_x):
        return mod_one(2 * raw_x) if self.is_bottom_side_moebius() else raw_x

    def adj_side(self, raw_point, raw_side):
        if not self.is_bottom_side_moebius():
            return raw_side
        return TOP if raw_point < 0.5 else BOTTOM

    def adj_to_raw(self, adj_x, adj_side):
        return adj_side * Rational('1/2') + adj_x/2 \
            if self.is_bottom_side_moebius() else adj_x






    




        
        
        

                

    def in_which_interval(self, point, side):
        r"""
        Finds the containing interval of a value on the specified side.

        INPUT:

        - ``value`` - a real number, in the interval $[0, 1)$

        - ``side`` - 0 or 1, 0 for the top side, 1 for the bottom. 

        OUTPUT:

        - non-negative integer -- the index of the containing interval
            on the specified side. If ``value`` is too close to one 
            of the division points, SaddleConnectionError is raised.
            Except when ``side`` is 1, and the bottom side is a 
            Moebius band, when it doesn't really matter whether 
            the point is in the 0th or 1st "not really existing"
            interval.

        TESTS::

            sage: f = Foliation.nonorientable_arnoux_yoccoz(4); f
            1 1 2 2 3 3
            Moebius band
            Lengths: (0.271844506346, 0.147798871261, 0.0803566223929)
            sage: f.in_which_interval(0.2, 0)
            0
            sage: f.in_which_interval(0.9, 0)
            4
            sage: f.in_which_interval(0.9999999999999, 0)
            Traceback (most recent call last):
            ...
            SaddleConnectionError
            sage: f.in_which_interval(0.271844506334, 0)
            Traceback (most recent call last):
            ...
            SaddleConnectionError
            sage: f.in_which_interval(0.499999999999, 1)
            0
            sage: f.in_which_interval(0.500000000001, 1)
            1

        """
        from bisect import bisect
        raw_point = self.adj_to_raw(point, side)
        raw_side = 0 if self.is_bottom_side_moebius() else side
        interval = bisect(self._divpoints[raw_side], raw_point)

        if interval == 0:
            to_check = {self._divpoints[raw_side][0]}
        elif interval == len(self._divpoints[raw_side]):
            to_check = {self._divpoints[raw_side][interval - 1],
                    self._divpoints[raw_side][0] + 1}
        else:
            to_check = {self._divpoints[raw_side][k]
                    for k in {interval - 1, interval}}

        if any(abs(raw_point - x) < EPSILON for x in to_check):
            # print self
            # print self._divpoints
            # print side, point
            # print raw_side, raw_point
            # print to_check
            raise SaddleConnectionError()
        return Interval(raw_side, (interval - 1) % 
                self.num_intervals(raw_side))


    def apply_iet(self, adj_point, adj_side, containing_interval):
        r"""
        Calculates the image of a point under the interval exchange map 
        and also the side and position of the image point.

        INPUT:

        - ``side`` - 0 or 1, the side the point is at (0 for the top,
            1 for the bottom side)

        - ``point`` - PointWithCoefficients, the point to map

        - ``containing_interval`` - the index of the interval containing
            the point on the given side. It may be omitted in which 
            case it is calculated. One reason for specifying it would be 
            saving time (if one knows it already). 

        OUTPUT:

            - tuple -- ((side, pos), point), where side is 0 or 1 and pos
                is the position of point which is a PointWithCoefficients

        TESTS::

            sage: f = Foliation.nonorientable_arnoux_yoccoz(4)
            sage: p = f._divpoints[1][1]; p
            (0.5, (0, 0, 0, 1, 1))
            sage: f.apply_iet(0, p, containing_interval = 1)
            ((0, 0), (0.228155493654, (-1, 0, 0, 1, 1)))
            sage: f.apply_iet(0, f._divpoints[0][0])
            Traceback (most recent call last):
            ...
            SaddleConnectionError

            sage: from sage.dynamics.foliations.foliation \
                    import PointWithCoefficients
            sage: f.apply_iet(0, PointWithCoefficients(0.9, [1,2,3,4,0]))
            ((0, 5), (0.980356622393, (1, 2, 4, 4, 0)))
            sage: f.apply_iet(1, PointWithCoefficients(0.9, [1,2,3,4,0]))
            ((1, 0), (0.4, (1, 2, 3, 3, 0)))

        """
        raw_point = self.adj_to_raw(adj_point, adj_side)
        new_int = containing_interval.pair(self)
        diff = raw_point - containing_interval.raw_endpoint(LEFT, self)
        if not containing_interval.is_flipped(self):
            raw_x = mod_one(new_int.raw_endpoint(LEFT, self) + diff)
        else:
            raw_x = mod_one(new_int.raw_endpoint(RIGHT, self) - diff)

        return (self.adj_side(raw_x, new_int.side), self.raw_to_adj(raw_x), new_int)

    # def point_int_on_other_side(self, point, side):
    #     if self.is_bottom_side_moebius():
    #         p = mod_one(point + Rational('1/2'))
    #         new_side = 0
    #     else:
    #         p = point
    #         new_side = (side + 1) % 2
    #     return (p, self.in_which_interval(p, new_side))

    def transform(self, interval, hdir):
        from separatrix import Separatrix
        from transition_map import new_foliation
        separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
        return new_foliation(separatrices, interval.endpoint(LEFT, self),
                             interval.endpoint_side(LEFT, self),
                             hdir = hdir,
                             is_one_sided = self.is_bottom_side_moebius(),
                             transformation_coding = SymmetryCoding(interval, hdir))
                             

    # def rotated(self, n):
    #     from separatrix import Separatrix
    #     from transition_map import new_foliation
    #     separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
    #     i = Interval(0,0).add_to_position(-n, self)
    #     return new_foliation(separatrices, i.endpoint(LEFT, self),
    #                          i.endpoint_side(LEFT, self),
    #                          is_one_sided = self.is_bottom_side_moebius())

    # def reversed(self):
    #     from separatrix import Separatrix
    #     from transition_map import new_foliation
    #     separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
    #     return new_foliation(separatrices, 0, TOP, direction = 'left',
    #                          is_one_sided = self.is_bottom_side_moebius())
                                 
    
    # def flip_over(self):
    #     if self.is_bottom_side_moebius():
    #         raise ValueError("Can't flip over a foliation around a "
    #                          "one-sided curve.")
    #     from separatrix import Separatrix
    #     from transition_map import new_foliation
    #     separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
    #     return new_foliation(separatrices, 0, BOTTOM, direction = 'right',
    #                          is_one_sided = False)
        

    def double_cover(self, foliation_or_surface):
        """Return an orienting double cover of the foliation.

        There are two types of orientation double covers. The first is
        the one that orients the surface, the second orients the
        foliation.

        Let's start with the orientation cover of the foliation.  If
        the foliation is already orientable (i.e. no interval is
        flipped), then there is nothing to do.

        So assume that it is not the case and there is at least one
        flipped pair of intervals. Our transverse curve has trivial
        holonomy, so it has two lifts. A leaf segment
        with endpoints in a not flipped pair of intervals lifts to 
        two segments, each having endpoints in one of the curves.
        A leaf segments with endpoints in a flipped pair of intervals
        lifts to two segments, one starting from the first curve, and
        ending in the second, the other the other way around.

        So a leaf segment with endpoints in the first lift
        of our curve is either just a lift of a segment corresponding
        to a not-flippped interval as above, or a concatenation of
        many segments, starting and ending with a segment that comes
        from a flipped interval, and with segments coming from
        not-flipped intervals in the middle. (One first travels from
        the first curve to the second, then returns a bunch of times
        to the second, and finally returns to the first.)

        So division points of the lifted foliation are the old 
        division points union the following points. From each 
        singularity travel along a separatrix, and stop just after
        going through a strip of a flipped pair of intervals.
        These endpoints are also new division points, so one has twice
        as many division points for the lift foliation as for the
        original one.
        
        Now turn to the double cover that orients the surface.
        If the surface is already orientable, there is nothing to do.
        If not, then our transverse curve, which separates the
        non-orientable surface to a Moebius band and another surface
        has two lifts. The Moebius band lifts to an annulus bounded
        by these two lifts. 
        
        If the other surface is orientable (i.e.
        all intervals are flipped),
        it has two lifts homeomorphic to itself, and they are
        glued to the two sides of the annulus. The foliation inside
        the annulus has the effect of a twist of 1/2 between the 
        surfaces on the two sides. All pairs of intervals remain
        flipped.

        If the complement of the Mobius band is non-orientable (i.e.
        there is at least one interval which is not flipped), then
        its double cover has one component with two bounding circles
        that are glued to the two boundary circles of the annulus.
        Arcs that were flipped stay flipped and will appear
        on both sides of the transverse curve. Arcs that were
        not flipped will turn into a pair of intervals on different
        sides, still not flipped (i.e. strips connecting 
        different sides of the annulus).

        In any case, one only has to change the interval pairings for
        the new involution, but the flips are inherited, and add a 
        1/2 twist.

        OUTPUT:

        - Foliation - if self is already orientable, the output is
          itself. Otherwise it is the foliation obtained by taking
          the double cover that orients the surface.

        EXAMPLES::
            
            sage: f = Foliation.nonorientable_arnoux_yoccoz(4)
            sage: g = Foliation.orientable_arnoux_yoccoz(3)
            sage: f.surface_orientable_double_cover() == g
            True
            
            sage: f = Foliation.nonorientable_arnoux_yoccoz(10)
            sage: g = Foliation.orientable_arnoux_yoccoz(9)
            sage: f.surface_orientable_double_cover() == g
            True

        """
        
        if self.is_orientable(foliation_or_surface):
            return self
        from separatrix import Separatrix
        from transition_map import new_foliation
        separatrices = Separatrix.get_all(self, number_of_flips_to_stop=0)
        if foliation_or_surface == 'foliation':
            separatrices += Separatrix.get_all(self, number_of_flips_to_stop=1)
            # removing removable singularities
            separatrices = [s for s in separatrices
                            if s.first_interval(0).prev(self) !=
                            s.first_interval(0).pair(self)]

        elif foliation_or_surface == 'surface':
            if not self.is_bottom_side_moebius():
                separatrices += Separatrix.get_all(self,
                                    stop_at_first_orientation_reverse=True)

        return new_foliation(separatrices, 0, 0,
                             is_one_sided = self.is_bottom_side_moebius(),
                             lift_type = foliation_or_surface)[0]
        

