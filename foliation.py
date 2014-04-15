from sage.dynamics.interval_exchanges.constructors import GeneralizedPermutation
from sage.structure.sage_object import SageObject
from collections import namedtuple
from sage.matrix.constructor import vector
from mymath import mod_one

from constants import *
from myexceptions import SaddleConnectionError


# def basis_vector(n, k):
#     """
#     Returns a standard basis vector of a vector space.

#     INPUT:

#     - ``n`` - positive integer, the dimension of the vector space

#     - ``k`` - 0 <= k < n, the index of the only coordinate which
#       is 1.

#     OUTPUT:

#     - list - the basis vector in as a list

#     EXAMPLES::

#         sage: from sage.dynamics.foliations.foliation import basis_vector
#         sage: basis_vector(5, 1)
#         [0, 1, 0, 0, 0]
#         sage: basis_vector(7, 6)
#         [0, 0, 0, 0, 0, 0, 1]
#     """
#     l = [0] * n
#     l[k] = 1
#     return l









from sage.rings.real_double import RDF







from sage.rings.rational import Rational
class Foliation(SageObject):
    """
    A measured foliation on a surface of finite type.

    Given a measured foliation we consider a two-sided simple closed
    curve on the surface which is transverse to the foliation. This 
    gives an interval exchange of the circle which is represented by
    an Involutions, length parameters, and a twist parameter. The 
    measure is always normalized so that the length of the curve is 1.

    We consider only surfaces with negative Euler characteristic, so
    foliations on the closed tori and closed Klein bottle are not 
    represented. By the Euler-Poincare theorem the foliation must have
    an even number of separatrices, say $2k$, where $k > 0$. There
    are two cases: all these separatrices hit our simple closed curve
    from above, or there is at least one separatrix hitting each side.

    In the first case the bottom side of the curve is a Moebius band,
    so we don't need a twist parameter. But even though we wouldn't
    need a length parameter for the bottom side to uniquely
    represent the foliation (it can be expressed in terms of the 
    length parameters above the curve), we consider 1/2 as the twist
    paratemeter for the bottom side with interval exchange imagined
    as 'z z'. Thus we have $k + 1$ parameters.

    In the other case there are $k$ length parameters (which again
    may not be independent), but a twist parameter here is 
    necessary which is $k + 1$ parameters again. This consistency may
    be useful for dealing with transition matrices between parameters
    of the same foliation, but looked at from the perspective of
    different simple closed curves.

    INPUT:

    - ``involution`` - an Involution, serving as the combinatorial
      description of the interval exchange

    - ``lenghts`` - as in iet.IntervalExchangeTransformation,
      this is either a list, tuple or dictionary of the length parameters.
      If it is a list or tuple, then the lengths are assigned to intervals
      in the order of their appearance in the Involution (their
      "index"). If it's a dict, then lengths should be assigned to
      the names of the intervals in the Involution.

    - ``twist`` - a real number, the twist parameter. This can 
      (and recommended to) be omitted if the bottom side is a 
      Moebius band, otherwise it is mandatory, because defaulting to
      zero would lead to an immediate saddle connection.

    EXAMPLES:
    
    Here are two different but equivalent definitions::

        sage: i = Involution('a a b b', 'c c', flips = 'b')
        sage: f = Foliation(i, [1, 2, 3], 1/2); f
        a a -b -b
        c c
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/12

        sage: g = Foliation(i, {'a':1,'b':2,'c':3}, 1/2); f
        a a -b -b
        c c
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/12

        sage: f == g
        True

    Omitting the twist if the bottom side is not a Moebius band
    throws an error::

        sage: Foliation(i, [1, 2, 3])
        Traceback (most recent call last):
        ...
        ValueError: The twist must be specified unless the bottom side is a Moebius band.

    When the bottom side is a Moebius band, it is okay to omit the 
    twist::

        sage: Foliation(Involution('a a b b c c'), [2, 2, 1])
        a a b b c c
        Moebius band
        Lengths: (1/5, 1/5, 1/10)

    The sum of the lengths on top and the sum at the bottom should
    (approximately) equal::

        sage: Foliation(i, [1, 1, 3], 1/2)
        Traceback (most recent call last):
        ...
        ValueError: The total length on the top and bottom are inconsistent.

    If they are just approximately equal, then on top of the
    normalization, the lengths are adjusted so the sums on top and 
    bottom are equal::

        sage: Foliation(i, [1, 2, 3.0000000000001], 1/2)
        a a -b -b
        c c
        Lengths: (0.166666666666661, 0.333333333333322, 0.500000000000017)
        Twist: 0.0833333333333306

    In reality, the same twisted interval exchange transformation can
    be represented by different Involutions and twist parameters. 
    For instance depending on which separatrix is chosen to be zero,
    the top letters can be 'a a b b', 'a b b a', 'b b a a' or
    'b a a b'. But even if one chooses one of these, varying the twist
    and the Involution can result in the same Foliation::

        sage: i = Involution('a b c', 'a c b'); i
        a b c
        a c b
        sage: f = Foliation(i, [1, 2, 3], 1); f
        a b c
        a c b
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/6

        sage: j = Involution('a b c', 'c b a'); j
        a b c
        c b a
        sage: g = Foliation(j, [1, 2, 3], 2); g
        a b c
        a c b
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/6

        sage: f == g
        True

    So rotate the bottom row of the Involution and change the twist
    if necessary to obtain the smallest possible positive twist. One
    can easily check that the Foliations f and g above have an 
    immediate saddle connection therefore they are not candidates for
    pseudo-anosov stable foliations. We don't check this in the
    constructor, only later when separatrices are lengthened to find
    other curves.

    """ 


    def __init__(self, top_letters, bottom_letters, lengths, flips = [], 
            twist = None):
        """
        TESTS::

            sage: i = Involution('a b c', 'a c b'); i
            a b c
            a c b
            sage: f = Foliation(i, [1, 2, 3], 1); f
            a b c
            a c b
            Lengths: (1/6, 1/3, 1/2)
            Twist: 1/6
            sage: f._divpoints
            [[(0, (0, 0, 0, 0)), (1/6, (1, 0, 0, 0)), (1/2, (1, 1, 0, 0))], [(1/6, (0, 0, 0, 1)), (1/3, (1, 0, 0, 1)), (5/6, (1, 0, 1, 1))]]
            sage: f._divvalues
            [[0, 1/6, 1/2], [1/6, 1/3, 5/6]]
            sage: f._lengths['a']
            (1/6, (1, 0, 0, 0))
            sage: f._lengths['b']
            (1/3, (0, 1, 0, 0))
            sage: f._lengths['c']
            (1/2, (0, 0, 1, 0))
            sage: f._twist
            (1/6, (0, 0, 0, 1))
            sage: f._involution
            a b c
            a c b

            sage: i = Involution('a a b b c c'); i
            a a b b c c
            Moebius band
            sage: f = Foliation(i, [1, 2, 3]); f
            a a b b c c
            Moebius band
            Lengths: (1/12, 1/6, 1/4)
            sage: f._divpoints
            [[(0, (0, 0, 0, 0, 0)),
              (1/12, (1, 0, 0, 0, 0)),
                (1/6, (2, 0, 0, 0, 0)),
                  (1/3, (2, 1, 0, 0, 0)),
                    (1/2, (2, 2, 0, 0, 0)),
                      (3/4, (2, 2, 1, 0, 0))],
                       [(0, (0, 0, 0, 0, 1)), (1/2, (0, 0, 0, 1, 1))]]
            sage: f._divvalues
            [[0, 1/12, 1/6, 1/3, 1/2, 3/4], [0, 1/2]]
            sage: f._lengths['a']
            (1/12, (1, 0, 0, 0, 0))
            sage: f._lengths['b']
            (1/6, (0, 1, 0, 0, 0))
            sage: f._lengths['c']
            (1/4, (0, 0, 1, 0, 0))
            sage: f._twist
            (0, (0, 0, 0, 0, 1))
            sage: f._involution
            a a b b c c
            Moebius band

        """
        self._tt = None # Train Track is created only when first used.

        if bottom_letters == 'moebius': # bottom side is Moebius
            bottom_letters = 'JOKER JOKER'
        self._init_gen_perm(top_letters, bottom_letters, flips)

        
        from bisect import bisect_left
        if not self.is_bottom_side_moebius():
            if twist == None:
                raise ValueError('The twist must be specified '
                'unless the bottom side is a Moebius band.')

        if isinstance(lengths, (list, tuple)):
            if len(lengths) != len(self.alphabet()):
                    raise ValueError('Bad number of lengths')
            self._lengths = {label: lengths[self._index_of_label[label]] 
                    for label in self.alphabet()}
                        
        if isinstance(lengths, dict):
            if set(self.alphabet()) != set(lengths.keys()):
                raise ValueError('Invalid length specification')     
            self._lengths = dict(lengths)

        if any(v <= 0 for v in self._lengths.values()):
            raise ValueError('Lengths must be positive')

        if self.is_bottom_side_moebius():
            self._lengths['JOKER'] = sum(self._lengths.values())
            twist = 0

        totals = [sum(interval.length() for interval in 
            self._all_intervals[side]) for side in {0,1}]

        if abs(totals[0] - totals[1]) > epsilon:
            raise ValueError('The total length on the top and '
                    'bottom are inconsistent.')
        
        #adjusting lengths in case they 
        #differ slightly on top/bottom
        for interval in self._all_intervals[0]:
            if interval.pair().side == interval.side:
                self._lengths[interval.label()] += \
                        (totals[1] - totals[0])/2
                break

        for label in self._lengths:
            self._lengths[label] /= totals[1]
        self._divvalues = [[0], [0]]

        for interval in self._all_intervals[0] + self._all_intervals[1]:
            self._divvalues[interval.side].append(self._divvalues[
                interval.side][-1] + self._lengths[interval.label()])
        for side in {0,1}:
            self._divvalues[side].pop()

        preimage_of_zero = mod_one(-twist/totals[1])
        containing_int = bisect_left(self._divvalues[1], 
                preimage_of_zero) % self.num_intervals(1)
        
        self._init_gen_perm(self._gen_perm_list[0],
                            self._gen_perm_list[1][containing_int:]+
                            self._gen_perm_list[1][:containing_int], 
                            flips = self.flips())
        self._twist = mod_one(self._divvalues[1][containing_int] -
                preimage_of_zero)
        self._divvalues[1] = [self._twist]
        for interval in self._all_intervals[1][:-1]:
            self._divvalues[1].append(self._divvalues[1][-1] + 
                    interval.length())
        
        self._length_twist_vector = [0] * len(self._lengths)
        self._length_twist_vector.append(self._twist)
        for label in self._lengths:
            self._length_twist_vector[self.index_of_label(label)]\
                    = self._lengths[label]
        self._length_twist_vector = vector(self._length_twist_vector)

        self._init_singularity_partition()

    def _init_singularity_partition(self):
        # initializing self._singularity_partition
        tt = self.train_track
        done = set()
        partition = []
        edge_circles = []
        for interval in self.intervals():
            if interval in done:
                continue
            new_interval = interval
            partition.append([])
            edge_circles.append([])
            end = 0 # we start from the left side of the interval
            while True:
                sing_index = new_interval if end == 0 else new_interval.next()
                if sing_index == interval and len(partition[-1]) > 0:
                    break
                done.add(sing_index)
                partition[-1].append(sing_index)
                edge_circles[-1].append(tt.get_center_edge(new_interval, end))
                end = (end + 1) % 2
                new_interval = new_interval.add_to_position((-1)**end)
                edge_circles[-1].append(tt.get_center_edge(new_interval, end).reversed())
                edge_circles[-1].append(tt.get_oriented_edge(new_interval,
                                            new_interval.pair(), 'pair'))
                new_interval = new_interval.pair()
                if new_interval.is_flipped():
                    end = (end + 1) % 2
                
        self._singularity_partition = partition
        self._paths_around_singularities = edge_circles
        
    @property
    def paths_around_singularities(self):
        return self._paths_around_singularities

    def _init_gen_perm(self, top_letters, bottom_letters, flips):
        self._gen_perm = GeneralizedPermutation(\
                        top_letters, bottom_letters, flips = flips)
        self._gen_perm_list = self._gen_perm.list()

        # initializing self._index_of_label and self._pair
        self._all_intervals = [[], []]
        label_to_interval = {}
        self._index_of_label = {}
        count = 0
        self._pair = {}
        for side in {0, 1}:
            for index in range(len(self._gen_perm_list[side])):
                interval = self.Interval(side, index, self)
                self._all_intervals[side].append(interval)
                label = interval.label()
                if label not in self._index_of_label:
                    self._index_of_label[label] = count
                    count += 1
                    label_to_interval[label] = interval
                else:
                    self._pair[label_to_interval[label]] = interval
                    self._pair[interval] = label_to_interval[label]


    def permutation(self):
        return self._gen_perm

    def intervals(self):
        if self.is_bottom_side_moebius():
            return self._all_intervals[0]
        return self._all_intervals[0] + self._all_intervals[1]

    def interval(self, side, index):
        return self._all_intervals[side][index]

    def num_intervals(self, side):
        return len(self._gen_perm[side])

    def labels(self):
        return self._gen_perm_list

    @property
    def train_track(self):
        from train_track import TrainTrack
        if not isinstance(self._tt, TrainTrack):
            self._tt = TrainTrack(self)
        return self._tt

    def index_of_label(self, label):
        """
        Returns the index of an letter.

        If n letters are used to notate the involution, they
        are indexed from 0 to n-1 according to their first 
        occurrence when read from to to bottom, from left to
        right.

        INPUT:

        - ``letter`` - string

        OUTPUT:

        - integer - the index of the letter

        EXAMPLES::

            sage: i = Involution('a b a c','d c d b')
            sage: i.index('a')
            0
            sage: i.index('b')
            1
            sage: i.index('c')
            2
            sage: i.index('d')
            3

        """
        return self._index_of_label[label]

    def raw_to_adj(self, raw_x):
        return mod_one(2 * raw_x) if self.is_bottom_side_moebius() else raw_x

    def _get_side(self, raw_point, raw_side):
        if not self.is_bottom_side_moebius():
            return raw_side
        return TOP if raw_point < 0.5 else BOTTOM

    def adj_to_raw(self, adj_x, adj_side):
        return adj_side * Rational('1/2') + adj_x/2 \
            if self.is_bottom_side_moebius() else adj_x

    class Interval(namedtuple('Interval', 'side, index')):
        def __new__(cls, side, index, foliation):
            self = super(Foliation.Interval, cls).__new__(cls, side, index)
            self._foliation = foliation
            return self
            
        def __repr__(self):
            return repr((self.side, self.index))
            
        def as_tuple(self):
            return (self.side, self.index)

        def add_to_position(self, n):
            return self._foliation._all_intervals[self.side][(self.index + n) 
                    % self._foliation.num_intervals(self.side)]

        def raw_endpoint(self, hdir):
            new_hdir = LEFT if hdir == MID else hdir
            x = self._foliation._divvalues[self.side][(self.index + new_hdir)
                    % self._foliation.num_intervals(self.side)]
            if hdir == MID:
                return mod_one(x + self.length()/2)
            else:
                return x

        def endpoint(self, hdir):
            return self._foliation.raw_to_adj(self.raw_endpoint(hdir))

        def endpoint_side(self, hdir):
            return self._foliation._get_side(self.raw_endpoint(hdir), self.side)

        def next(self):
            return self.add_to_position(1)

        def prev(self):
            return self.add_to_position(-1)

        def label(self):
            return self._foliation.labels()[self.side][self.index]

        def length(self):
            return self._foliation._lengths[self.label()]

        def is_wrapping(self):
            return self.endpoint(LEFT) > self.endpoint(RIGHT)

        def is_orienation_reversing(self):
            return self.is_flipped() != (self.pair().side == self.side)

        def is_flipped(self):
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
            x = self._foliation._gen_perm[self.side][self.index]
            if isinstance(x, tuple):
                # self._gen_perm is a FlippedLabelledPermutationLI
                return x[1] == -1
            #self._gen_perm is a LabelledPermutationIET 
            return False


        def pair(self):
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
            side, index = self._foliation._pair[(self.side, self.index)]
            return self._foliation._all_intervals[side][index]

        def which_singularity(self):
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
            
v            OUTPUT:

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
            sp = self._foliation._singularity_partition
            for i in range(len(sp)):
                if (self.side, self.index) in sp[i]:
                    return i
            raise ValueError("Invalid singularity specification.")





    def is_bottom_side_moebius(self):
        """
        Decides if the bottom side is Moebius band without
        punctures.

        This happends exactly when the second argument in
        the constructor is omitted.

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: i = Involution('a a b b')
            sage: i.is_bottom_side_moebius()
            True

        Here the bottom component is a once punctured Moebius
        band, so it returns False::

            sage: i = Involution('a a b b', 'c c')
            sage: i.is_bottom_side_moebius()
            False

        """
        return self._all_intervals[1][0].label() == 'JOKER'


    def flips(self):
        """
        Returns the list of flips.

        OUTPUT:

        - list -- the list of flipped interval names

        EXMAPLES::

            sage: i = Involution('a a b b','c c',flips = 'ab')
            sage: i.flips()
            ['a', 'b']

            sage: i = Involution('a a', flips = 'a')
            sage: i.flips()
            ['a']

        """
        if isinstance(self._gen_perm[0][0], tuple):
            # self._gen_perm is a FlippedLabelledPermutationLI
            # that has a flips() method
            return self._gen_perm.flips()
        else: #self._gen_perm is a LabelledPermutationIET 
            # hence there are no flips
            return []

    def alphabet(self):
        """
        Returns the set of interval names.

        OUTPUT:

        set -- the set of interval names

        EXAMPLES::

            sage: i = Involution('a a b b','c c',flips='c')
            sage: i.alphabet() == {'a', 'b', 'c'}
            True

        """
        s = set(self._gen_perm.alphabet())
        s.discard('JOKER')
        return s

    def is_orientable(self, foliation_or_surface):
        """
        Decides if the suspension foliations are orientable.

        OUTPUT:

        - boolean --

        EXAMPLES:

        If there are flipped intervals, the foliation is 
        non-orientable::

            sage: i = Involution('a a b b', 'c c', flips ='a')
            sage: i.is_foliation_orientable()
            False

        If there are no flips, the foliation is orientable::

            sage: i = Involution('a a b b', 'c c')
            sage: i.is_foliation_orientable()
            True

        Decides if the suspension surface is orientable.

        OUTPUT:

        - boolean -- 

        EXAMPLES:

        If the bottom side is a Moebius band, it is always
        non-orientable::

            sage: i = Involution('a a b b')
            sage: i.is_surface_orientable()
            False

        Or if there is flipped pair on different sides::

            sage: i = Involution('a b c','c b a', flips = 'a')
            sage: i.is_surface_orientable()
            False

        Or if there is not flipped pair on the same side::

            sage: i = Involution('a a b b', 'c c', flips='ac')
            sage: i.is_surface_orientable()
            False

        Otherwise it is orientable::

            sage: i = Involution('a b c', 'b c a')
            sage: i.is_surface_orientable()
            True

        """
        if foliation_or_surface == 'foliation':
            return len(self.flips()) == 0
        if foliation_or_surface == 'surface':
            return all(not i.is_orienation_reversing() for side in range(2)
                       for i in self._all_intervals[side])
            # including the Moebius side as well


    def singularity_partition(self):
        """
        Returns the singularity partition of self.

        OUTPUT:

        - list of lists - each element of the list is a list corresponding to
          a singularity, and each such list contains the tuples of positions
          that are being identified

        EXAMPLES::

            sage: i = Involution('a a b b c c')
            sage: sp = i.singularity_partition()
            sage: len(sp) == 1
            True
            sage: set(sp[0]) == {(0,2), (0,3), (0, 4), (0, 5), (0, 0), (0, 1)}
            True

            sage: i = Involution('a b c d', 'b a d c', flips = 'ac')
            sage: sp = i.singularity_partition()
            sage: len(sp) == 2
            True
            sage: set(sp[0]) == {(1, 1), (0, 2), (1, 0), (0, 1)}
            True
            sage: set(sp[1]) == {(0, 0), (1, 3), (0, 3), (1, 2)}
            True

        """
        return list(self._singularity_partition)

         
    def singularity_type(self):
        """
        Returns the singularity type of self.

        The suspension of the Involution yields a foliation.
        The singularity type of that foliation is the tuple of
        the number of prongs at singularities.

        OUTPUT:

        - tuple -

        EXAMPLES::

            sage: i = Involution('a a b b c c', flips = 'abc'); i
            -a -a -b -b -c -c
            Moebius band
            sage: i.singularity_type()
            (3, 1, 1, 1)

            sage: i = Involution('a a b b c c', 'd d', flips = 'abc'); i
            -a -a -b -b -c -c
             d  d
            sage: i.singularity_type()
            (3, 2, 1, 1, 1)

            sage: i = Involution('a', 'a'); i
            a
            a
            sage: i.singularity_type()
            (2,)

            sage: i = Involution('a b', 'a b'); i
            a b 
            a b
            sage: i.singularity_type()
            (2, 2)

        """
        t = sorted([x for x in map(len, 
            self._singularity_partition)], reverse = True)
        return tuple(t)

    @property
    def divvalues(self):
        return self._divvalues

        

    def with_changed_lengths(self, length_vector):
        if self.is_bottom_side_moebius():
            return Foliation(self._gen_perm_list[0], 'moebius', length_vector,
                             flips = self.flips())
        return Foliation(self._gen_perm_list[0], self._gen_perm_list[1],
                         length_vector[:-1], flips = self.flips(),
                         twist = length_vector[-1])

    def __eq__(self, other):
        r"""
        Tests approximate equality.

        The involutions have to equal up to different labelling, but the
        length and twist parameters are enough to be very close for the
        equality of Foliations.

        EXAMPLES::

            sage: i = Involution('a a b b', 'c c')
            sage: f = Foliation(i, [1, 2, 3], 1/2); f
            a a b b
            c c
            Lengths: (1/6, 1/3, 1/2)
            Twist: 1/12

            sage: j = Involution('1 1 2 2', '3 3')
            sage: g = Foliation(j, [1, 2, 3.000000000001], 1/2); g
            1 1 2 2
            3 3
            Lengths: (0.166666666666611, 0.333333333333222, 0.500000000000167)
            Twist: 0.0833333333333055

            sage: f == g
            True

        """
        return self.equals(other, epsilon)

    def equals(self, other, allowed_error):
        return self.is_bottom_side_moebius() == other.is_bottom_side_moebius()\
                and self._gen_perm == other._gen_perm and \
                abs(self._length_twist_vector -
                    other._length_twist_vector) < allowed_error
        

    def _repr_(self):
        r"""
        Returns a representation of self.

        TESTS::

            sage: Foliation(Involution('a b c', 'a c b'), [1, 2, 3], 1)
            a b c
            a c b
            Lengths: (1/6, 1/3, 1/2)
            Twist: 1/6

            sage: Foliation(Involution('a a b b c c'), [2, 2, 1])
            a a b b c c
            Moebius band
            Lengths: (1/5, 1/5, 1/10)

        """
        if self.is_bottom_side_moebius():
            d = dict(self._lengths)
            del d['JOKER']
            return "{0}\nLengths: {1}".format(\
                repr(self._gen_perm).split('\n')[0] + '\nMoebius band', d)
        return "{0}\nLengths: {1}\nTwist: {2}".format(repr(self._gen_perm),
                self._lengths, self._twist)

    def _latex_(self):
        return self.latex_options().tikz_picture()

    def latex_options(self):
        if not hasattr(self, '_latex_opts'):
            from foliation_latex import FoliationLatex
            self._latex_opts = FoliationLatex(self)
        return self._latex_opts

    def set_latex_options(self, **kwds):
        opts = self.latex_options()
        opts.set_options(**kwds)

                

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
        interval = bisect(self._divvalues[raw_side], raw_point)

        if interval == 0:
            to_check = {self._divvalues[raw_side][0]}
        elif interval == len(self._divvalues[raw_side]):
            to_check = {self._divvalues[raw_side][interval - 1],
                    self._divvalues[raw_side][0] + 1}
        else:
            to_check = {self._divvalues[raw_side][k]
                    for k in {interval - 1, interval}}

        if any(abs(raw_point - x) < epsilon for x in to_check):
            # print self
            # print self._divvalues
            # print side, point
            # print raw_side, raw_point
            # print to_check
            raise SaddleConnectionError()
        return self._all_intervals[raw_side][(interval - 1) % 
                self.num_intervals(raw_side)]


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
        new_int = containing_interval.pair()
        diff = raw_point - containing_interval.raw_endpoint(LEFT)
        if not containing_interval.is_flipped():
            raw_x = mod_one(new_int.raw_endpoint(LEFT) + diff)
        else:
            raw_x = mod_one(new_int.raw_endpoint(RIGHT) - diff)

        return (self._get_side(raw_x, new_int.side), self.raw_to_adj(raw_x), new_int)

    # def point_int_on_other_side(self, point, side):
    #     if self.is_bottom_side_moebius():
    #         p = mod_one(point + Rational('1/2'))
    #         new_side = 0
    #     else:
    #         p = point
    #         new_side = (side + 1) % 2
    #     return (p, self.in_which_interval(p, new_side))


    def rotated(self, n):
        from separatrix import Separatrix
        from transition_map import new_foliation, get_tt_map
        separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
        i = self._all_intervals[0][0]
        i = i.add_to_position(-n)
        return new_foliation(separatrices, i.endpoint(LEFT), i.endpoint_side(LEFT),
                             is_one_sided = self.is_bottom_side_moebius())

    def reversed(self):
        from separatrix import Separatrix
        from transition_map import new_foliation, get_tt_map
        separatrices = Separatrix.get_all(self, number_of_flips_to_stop = 0)
        return new_foliation(separatrices, 0, TOP, direction = 'left',
                             is_one_sided = self.is_bottom_side_moebius())
                                 
    


    def double_cover(self, foliation_or_surface):
        """
        Returns the orienting double cover of the foliation.

        If the foliation is already orientable (i.e. the surface
        is orientable and the holonomy is trivial, or the surface
        is non-orientable and the holonomy is Z_2, or equivalently
        no interval is flipped), then there is nothing to do.

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

        Returns the double cover of the foliation that orients the
        surface.

        If the surface is already orientable, there is nothing to do.
        If not, then the our transverse curve, which separates the
        non-orientable surface to a Moebius band and another surface
        has two lifts. The Moebius band lifts to an annulus bounded
        by these two lifts. 
        
        If the other surface is orientable (i.e.
        all intervals are flipped),
        it has two lifts homeomorphic to itself, and they are
        glued to the two sides of the annulus. The folition inside
        the annulus has the effect of a twist of 1/2 between the 
        surfaces on the two sides. All pairs of intervals remain
        flipped.

        If the complement of the Mobius band is non-orientable (i.e.
        there is at least one interval which is not twisted), then
        its double cover has one components with two bounding circles
        that are glued to the two boundary circles of the annulus.
        Arcs that were flipped stay flipped and will appear
        on both sides of the transverse curve. Arcs that were
        not flipped will turn into a pair of intervals on different
        sided, still not flipped (i.e. strips connecting 
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
                            if s.first_interval(0).prev() !=
                            s.first_interval(0).pair()]

        elif foliation_or_surface == 'surface':
            if not self.is_bottom_side_moebius():
                separatrices += Separatrix.get_all(self,
                                    stop_at_first_orientation_reverse=True)

        return new_foliation(separatrices, 0, 0,
                             is_one_sided = self.is_bottom_side_moebius(),
                             lift_type = foliation_or_surface)[0]
        

