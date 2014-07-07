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

from base import mod_one, LEFT, MID, RIGHT
from collections import namedtuple


class Interval(namedtuple('Interval', 'side, index')):

    def __repr__(self):
        return repr((self.side, self.index))

    def as_tuple(self):
        return (self.side, self.index)

    def add_to_position(self, n, fol):
        return Interval(self.side, (self.index + n) 
                        % fol.num_intervals(self.side))

    def raw_endpoint(self, hdir, fol):
        new_hdir = LEFT if hdir == MID else hdir
        x = fol._divpoints[self.side][(self.index + new_hdir)
                % fol.num_intervals(self.side)]
        if hdir == MID:
            return mod_one(x + self.length(fol)/2)
        else:
            return x

    def endpoint(self, hdir, fol):
        return fol.raw_to_adj(self.raw_endpoint(hdir, fol))

    def endpoint_side(self, hdir, fol):
        return fol.adj_side(self.raw_endpoint(hdir, fol), self.side)

    def next(self, fol):
        return self.add_to_position(1, fol)

    def prev(self, fol):
        return self.add_to_position(-1, fol)

    def label(self, fol):
        return fol._gen_perm_list[self.side][self.index]

    def numerical_label(self, fol):
        """Return the numerical label of the letter.

        The `n` labels used to code the permutation of intervals are
        associated with numerical labels from 0 to `n-1` to make some
        calculations more convenient. The numbers from 0 to `n-1` are
        associated to the labels in the order of their appearance.

        INPUT:

        - ``fol`` -- the underlying ``Foliation``

        OUTPUT:

        a non-negative integer, the numerical label of the interval

        EXAMPLES::

        sage: f = Foliation('a b b c a c','d d',flips='abc')
        sage: [Interval(0,i).numerical_label(f) for i in range(6)]
        [0, 1, 1, 2, 0, 2]
        sage: [Interval(1,i).numerical_label(f) for i in range(2)]        
        [3, 3]

        """
        return fol._numerical_label[self.label(fol)]
    
    def length(self, fol):
        return fol._lengths[self.label(fol)]

    def is_wrapping(self, fol):
        return self.endpoint(LEFT, fol) > self.endpoint(RIGHT, fol)

    def is_orientation_reversing(self, fol):
        return self.is_flipped(fol) != (self.pair(fol).side == self.side)

    def is_flipped(self, fol):
        """
        Decides if the interval at a certain position is 
        flipped.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.

        OUTPUT:

        - boolean -- True is the interval is flipped, False
          is not

        EXAMPLES::

            sage: i = Involution('a a b b','c c', flips='bc');i
            a a -b -b
            -c -c
            sage: i.is_flipped((0,0))
            False
            sage: i.is_flipped((0,1))
            False
            sage: i.is_flipped((0,2))
            True
            sage: i.is_flipped((0,3))
            True
            sage: i.is_flipped((1,0))
            True
            sage: i.is_flipped((1,1))
            True

        """
        x = fol._gen_perm[self.side][self.index]
        if isinstance(x, tuple):
            # self._gen_perm is a FlippedLabelledPermutationLI
            return x[1] == -1
        #self._gen_perm is a LabelledPermutationIET 
        return False


    def pair(self, fol):
        """
        Returns the position of the pair of the interval at
        a specified position.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.

        OUTPUT:

        - tuple -- the position of the pair

        EXAMPLES::

            sage: i = Involution('a b a b','c c', flips = 'ab')
            sage: i.pair((0,0))
            (0, 2)
            sage: i.pair((0,2))
            (0, 0)
            sage: i.pair((1,1))
            (1, 0)

        """
        side, index = fol._pair[(self.side, self.index)]
        return Interval(side, index)

    def which_singularity(self, fol):
        """
        Returns the index of the singularity for the beginning of each
        interval.

        There is a singularity of the foliation on the leaf containing the
        left endpoint of each interval for any suspension. There may be only
        one singularity of all the vertices are identified, or more if not.
        If there are $n$ singularities ($n\ge 1$), we assign 0, 1, ..., $n-1$
        to them in some order, this is called its index. The index is 
        therefore well-defined only up to a permutation of these values.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.

        OUTPUT:

        - integer - the index of the specified singularity. 

        EXAMPLES:

        The following Involution has 2 singularities, one has 5 prongs, the
        other 1 prong. The position of the 1-prong singularity is at (1,0).
        Here is a possible output:

            sage: i = Involution('a a b', 'c b c', flips = 'c')
            sage: i.singularity_type()
            (5, 1)
            sage: i.which_singularity((0,0)) 
            0
            sage: i.which_singularity((0,1)) 
            0
            sage: i.which_singularity((0,2)) 
            0
            sage: i.which_singularity((1,0)) 
            1
            sage: i.which_singularity((1,1)) 
            0
            sage: i.which_singularity((1,2)) 
            0

        """
        sp = fol.singularity_partition()
        for i in range(len(sp)):
            if (self.side, self.index) in sp[i]:
                return i
        raise ValueError("Invalid singularity specification.")


