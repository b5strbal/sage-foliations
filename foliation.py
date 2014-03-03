from sage.dynamics.interval_exchanges.constructors import GeneralizedPermutation
from sage.structure.sage_object import SageObject
from collections import namedtuple
from sage.matrix.constructor import vector
epsilon = 1e-10


def mod_one(x):
    """
    Returns a number modulo 1.

    INPUT:

    - ``x`` - a real number

    OUTPUT:

    - a real number of the same type as the input

    TESTS::

        sage: from sage.dynamics.foliations.foliation import mod_one
        sage: mod_one(2.5)
        0.500000000000000
        sage: mod_one(-1.7)
        0.300000000000000
        sage: mod_one(7/6)
        1/6
        sage: mod_one(-1/6)
        5/6
        sage: a = QQbar(sqrt(2)); a
        1.414213562373095?
        sage: mod_one(a)
        0.4142135623730951?

    """
    from sage.functions.other import floor 
    return x - floor(x)

def basis_vector(n, k):
    """
    Returns a standard basis vector of a vector space.

    INPUT:

    - ``n`` - positive integer, the dimension of the vector space

    - ``k`` - 0 <= k < n, the index of the only coordinate which
      is 1.

    OUTPUT:

    - list - the basis vector in as a list

    EXAMPLES::

        sage: from sage.dynamics.foliations.foliation import basis_vector
        sage: basis_vector(5, 1)
        [0, 1, 0, 0, 0]
        sage: basis_vector(7, 6)
        [0, 0, 0, 0, 0, 0, 1]

    """
    l = [0] * n
    l[k] = 1
    return l

def make_positive(vec):
    """
    Returns a parallel vector to a given vector with all positive
    coordinates if possible.

    INPUT:

    - ``vec`` - a vector or list or tuple of numbers

    OUTPUT:

    - vector - a vector with positive coordinates if this is possible,
      otherwise -1

    EXAMPLES:

    If all coordinates of the input vector are already positive, the
    same vector is returned::

        sage: from sage.dynamics.foliations.foliation import make_positive
        sage: v = vector([1, 3, 5])
        sage: make_positive(v) == v
        True

    If all are negative, its negative is returned::

        sage: make_positive((-1, -5))
        (1, 5)

    Even if the coordinates are complex, but have very small imaginary
    part as a result of an approximate eigenvector calculation for 
    example, the coordinates are treated as real::

        sage: make_positive((40.24 - 5.64e-16*I, 1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)
        sage: make_positive((-40.24 - 5.64e-16*I, -1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)

    If there is a complex coordinate which is not negligible, -1 is
    returned::

        sage: make_positive((-40.24 - 5.64e-6*I, -1.2 + 4.3e-14*I))
        -1

    If one coordinate is zero, or very close to zero, -1 is returned::

        sage: make_positive((1, 0, 2))
        -1
        sage: make_positive((-40.24e-15*I - 5.64e-16*I, -1.2))
        -1

    If there are both negative and positive coordinates, -1 is 
    returned::

        sage: make_positive((-3, 4, 5))
        -1
        
    """
    if any(abs(x) < epsilon for x in vec) or \
        any(abs(x.imag()) > epsilon for x in vec):
            return -1
    newvec = vector([x.real() for x in vec])
    if vec[0].real() < 0:
        newvec = -newvec
    if any(x < 0 for x in newvec):
        return -1
    return newvec

def get_good_eigendata(transition_matrix, is_twisted):
    """
    
        sage: from sage.dynamics.foliations.foliation import get_good_eigendata


    """
    ev = transition_matrix.eigenvectors_right()
    ret_list = []
    for x in ev:
        if abs(x[0].imag()) < epsilon and x[0].real() > 0 \
                and abs(x[0].real() - 1) > epsilon:
            for vec in x[1]:
                newvec = make_positive(vec)
                if newvec != -1:
                    norm = sum([abs(y) for y in newvec])
                    if is_twisted:
                        norm -= abs(newvec[-1])
                    normalized_vec = [abs(y / norm) for y in newvec]
                    ret_list.append((x[0].real(), normalized_vec))
    return ret_list













#class PointWithCoefficients(namedtuple("Pwc", 
#        "value, coefficients")):
#    """
#    Represents a real number together with its coefficients as
#    a linear combination of the a certain set of real numbers.
#
#    In practice this set is the set of interval lengths and 
#    possibly the twist of a Foliation.
#
#    INPUT:
#
#    - ``value`` - a real number (can be floating point, 
#      or an exact algebraic number)
#
#    - ``coefficients`` - a list or tuple or vector of 
#      the coefficients
#
#    EXAMPLES:
#
#    The coefficients can be specified as a list, tuple or 
#    vector::
#
#        sage: from sage.dynamics.foliations.foliation import \
#                PointWithCoefficients
#        sage: a = PointWithCoefficients(1.2, (3, -1, 0))
#        sage: b = PointWithCoefficients(0.8, [2, 1, 5])
#        sage: c = PointWithCoefficients(3.4, vector((2, 3)))
#
#    One can add or subtract two objects as long as they are 
#    of the same type:
#
#        sage: a + b
#        (2.00000000000000, (5, 0, 5))
#        sage: a - b
#        (0.400000000000000, (1, -2, -5))
#    
#    """
#    def __new__(cls, value, coefficients):
#        self = super(PointWithCoefficients, cls).__new__(cls, 
#                value, vector(coefficients))
#        return self
#
#    def __repr__(self):
#        """
#        Returns the representation of self.
#
#        TESTS::
#
#            sage: from sage.dynamics.foliations.foliation import \
#                PointWithCoefficients
#            sage: PointWithCoefficients(3, (4, 3, 2))
#            (3, (4, 3, 2))
#
#        """
#        return repr((self.value, self.coefficients))
#
#    def __add__(self, other):
#        """
#        Adds the numbers and their coefficient vectors.
#
#        TESTS::
#
#            sage: from sage.dynamics.foliations.foliation import PointWithCoefficients
#            sage: a = PointWithCoefficients(1.2, (3, -1, 0))
#            sage: b = PointWithCoefficients(0.8, (2, 1, 5))
#            sage: a + b
#            (2.00000000000000, (5, 0, 5))
#
#        """
#        return PointWithCoefficients(self.value + other.value,
#                self.coefficients + other.coefficients)
#        
#    def __sub__(self, other):
#        """
#        Subtracts the numbers and their coefficient vectors.
#
#        TESTS::
#
#            sage: from sage.dynamics.foliations.foliation import PointWithCoefficients
#            sage: a = PointWithCoefficients(1.2, (3, -1, 0))
#            sage: b = PointWithCoefficients(0.8, [2, 1, 5])
#            sage: a - b
#            (0.400000000000000, (1, -2, -5))
#        """
#        return PointWithCoefficients(self.value - other.value,
#                self.coefficients - other.coefficients)
#
#
#    def mod_one(self):
#        """
#        Returns the PointWithCoefficients corresponding to the
#        real number of self modulo 1.
#
#        The sum of the numbers in the generating set used for
#        linear combinations is 1 hence all but the twist
#        coefficient is decreased by the floor of the real
#        number of self.
#
#        OUTPUT:
#
#        - PointWithCoefficients --
#
#        EXAMPLES::
#
#            sage: from sage.dynamics.foliations.foliation import PointWithCoefficients
#            sage: p = PointWithCoefficients(3.2, (4, 5, 3))
#            sage: p.mod_one()
#            (0.200000000000000, (1, 2, 3))
#
#            sage: q = PointWithCoefficients(-1.3, (0, 2, 4,-2))
#            sage: q.mod_one()
#            (0.700000000000000, (2, 4, 6, -2))
#
#        """
#        from sage.functions.other import floor 
#        if self.value == 0: #need to check beacuse of a bug
#            # with algebraic numbers
#            n = 0
#        else:
#            n = floor(self.value)
#        return PointWithCoefficients(self.value - n, [x - n 
#            for x in self.coefficients[:-1]] + \
#                    [self.coefficients[-1]])



from sage.rings.real_double import RDF
def arnoux_yoccoz_factor(genus, field = RDF):
    """
    Returns the Perron root of the polynomials 
    $x^g - x^{g-1} - ... - 1$.

    INPUT:

    - ``genus`` - integer, the genus of the orientable
      closed Arnoux-Yoccoz surface, which is the same as the degree
      of the minimalpolynomial.

    - ``field`` - the number field in which roots of the polynomial
      are searched. The default is RDF, but one want the exact number
      in QQbar.

    OUTPUT:

    - number -- number type depends on the field specified

    EXAMPLES::

        sage: from sage.dynamics.foliations.foliation import arnoux_yoccoz_factor
        sage: arnoux_yoccoz_factor(3)
        1.83928675521
        sage: arnoux_yoccoz_factor(3, field = QQbar)
        1.839286755214161?

    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import \
        PolynomialRing
    R = PolynomialRing(ZZ, 't')
    poly = R([-1] * genus + [1])
    return max([abs(x[0]) for x in poly.roots(field)])








class SaddleConnectionError(Exception):
    pass

class RestrictionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value






    











 
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

        if bottom_letters == 'moebius': # bottom side is Moebius
            bottom_letters = 'JOKER JOKER'
        self._gen_perm = GeneralizedPermutation(\
                top_letters, bottom_letters, flips = flips)

        # initializing self._index_of_label and self._pair
        self._all_intervals = [[], []]
        label_to_interval = {}
        self._index_of_label = {}
        count = 0
        self._pair = {}
        for side in {0, 1}:
            for index in range(len(self.labels()[side])):
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


        # initializing self._singularity_partition
        done = set()
        partition = []
        for interval in self.intervals():
            if interval in done:
                continue
            new_interval = interval
            partition.append([])
            direction = 'left'
            while True:
                if direction == 'left':
                    new_interval = new_interval.prev().pair()
                    if not new_interval.is_flipped():
                        new_interval = new_interval.next()
                        direction = 'away'
                else:
                    new_interval = new_interval.pair()
                    if not new_interval.is_flipped():
                        direction = 'left'
                    else:
                        new_interval = new_interval.next()
                partition[-1].append(new_interval)
                done.add(new_interval)
                if interval == new_interval:
                    break

        self._singularity_partition = partition




        
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
            #self._half = PointWithCoefficients(Rational('1/2'),
            #        basis_vector(len(lcopy) + 1, len(lcopy) - 1))

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
        #for letter in lcopy:
        #    self._lengths[letter] = PointWithCoefficients(\
        #            lcopy[letter]/totals[1], 
        #            basis_vector(len(lcopy) + 1,
        #                involution.index(letter)))

        #self._divpoints = [[PointWithCoefficients(0,
        #    [0] * (len(lcopy) + 1))]
        #    for i in range(2)] 
        self._divvalues = [[0], [0]]

        #for i in range(2):
        #    for j in range(len(involution[i]) - 1):
        #        self._divpoints[i].append(self._divpoints[i][-1] +
        #                self._lengths[involution[i][j]])
        for interval in self._all_intervals[0] + self._all_intervals[1]:
            self._divvalues[interval.side].append(self._divvalues[
                interval.side][-1] + self._lengths[interval.label()])
        for side in {0,1}:
            self._divvalues[side].pop()

        #self._divvalues = [[x.value for x in self._divpoints[i]] 
        #    for i in range(2)]

        preimage_of_zero = mod_one(-twist/totals[1])
        containing_int = bisect_left(self._divvalues[1], 
                preimage_of_zero) % self.num_intervals(1)
        self._gen_perm = self._rotated_gen_perm(0, -containing_int)
        #self._twist = PointWithCoefficients(mod_one(\
        #        self._divpoints[1][containing_int].value - 
        #        preimage_of_zero), [0] * len(lcopy) + [1])
        self._twist = mod_one(self._divvalues[1][containing_int] -
                preimage_of_zero)
        #self._divpoints[1] = [self._twist]
        self._divvalues[1] = [self._twist]
        #for j in range(len(involution[1]) - 1):
        #    self._divpoints[1].append(self._divpoints[1][-1] +
        #            self._lengths[self._involution[1][j]])

        #self._divvalues[1] = [x.value for x in self._divpoints[1]] 
        for interval in self._all_intervals[1]:
            self._divvalues[1].append(self._divvalues[1][-1] + 
                    interval.length())

        self._length_twist_vector = [0] * len(self._lengths)
        self._length_twist_vector.append(self._twist)
        for label in self._lengths:
            self._length_twist_vector[self.index_of_label(label)]\
                    = self._lengths[label]
        self._length_twist_vector = vector(self._length_twist_vector)

    def intervals(self):
        if self.is_bottom_side_moebius():
            return self._all_intervals[0]
        return self._all_intervals[0] + self._all_intervals[1]

    def num_intervals(self, side):
        return len(self._gen_perm[side])

    def labels(self):
        return self._gen_perm.list()

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

    class Interval(namedtuple('Interval', 'side, index')):
        def __new__(cls, side, index, foliation):
            self = super(Foliation.Interval, cls).__new__(cls, side, index)
            self._foliation = foliation
            return self

        def add_to_position(self, n):
            return self._foliation._all_intervals[self.side][(self.index + n) 
                    % self._foliation.num_intervals(self.side)]

        def endpoint(self, pos):
            return self._foliation._divvalues[self.side][(self.index + pos)
                    % self._foliation.num_intervals(self.side)]

        def midpoint(self):
            return mod_one(self.endpoint(0) + self.length()/2)

        def next(self):
            return self.add_to_position(1)

        def prev(self):
            return self.add_to_position(-1)

        def label(self):
            return self._foliation.labels()[self.side][self.index]

        def length(self):
            return self._foliation._lengths[self.label()]

        def is_wrapping(self):
            return self.endpoint(0) > self.endpoint(1)


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
            return Foliation.Interval(side, index, self._foliation)

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

    def is_foliation_orientable(self):
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

        """
        return len(self.flips()) == 0

    def is_surface_orientable(self):
        """
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
        for interval in flatten(self._all_intervals):
            cond1 = (interval.side == interval.pair().side)
            cond2 = interval.is_flipped()
            if cond1 != cond2:
                return False
        return True

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




    @classmethod
    def orientable_arnoux_yoccoz(self, genus):
        r"""
        Constructs an orientable Arnoux-Yoccoz foliation.

        INPUT:

        - ``genus`` - integer, at least 3, the genus of the surface

        OUTPUT:

        - Foliation -- the Arnoux-Yoccoz foliation of specified genus

        EXAMPLES::

            sage: Foliation.orientable_arnoux_yoccoz(3)
            1 2 3 4 5 6
            4 3 6 5 2 1
            Lengths: (0.271844506346, 0.271844506346, 0.147798871261, 
                    0.147798871261, 0.0803566223929, 0.0803566223929)
            Twist: 0.0436890126921

        """
        if genus < 3:
            raise ValueError('The genus of an orientable'
                    'Arnoux-Yoccoz surface is at least 3')
        n = 2 * genus
        top = range(1, n + 1)
        def switch(k):
            if k % 2 == 0:
                return k + 1
            return k - 1
        bottom = [top[switch(i)] for i in range(n)]
        sf = arnoux_yoccoz_factor(genus)
        l = [1 / sf**(i + 1) for i in range(genus)] * 2
        return Foliation(top, bottom, sorted(l, reverse = True), twist = 1)

    @classmethod
    def nonorientable_arnoux_yoccoz(self, genus):
        r"""
        Constructs an Arnoux-Yoccoz foliation on a non-orientable surface.

        INPUT:

        - ``genus`` - integer, at least 4, the genus of the surface

        OUTPUT:

        - Foliation -- the Arnoux-Yoccoz foliation of specified genus

        EXAMPLES::

            sage: Foliation.nonorientable_arnoux_yoccoz(4)
            1 1 2 2 3 3
            Moebius band
            Lengths: (0.271844506346, 0.147798871261, 0.0803566223929)

        """
        if genus < 4:
            raise ValueError('The genus of a non-orientable '
                    'Arnoux-Yoccoz surface is at least 4')
        top = sorted(2 * range(1, genus))
        sf = arnoux_yoccoz_factor(genus - 1)
        return Foliation(top, 'moebius',[1/sf**i for i in range(genus - 1)])

    @classmethod
    def RP2_arnoux_yoccoz(self):
        r"""
        Constructs the Arnoux-Yoccoz foliation on RP2.

        OUTPUT:

        - Foliation -- the Arnoux-Yoccoz foliation on RP2

        EXAMPLES::

            sage: Foliation.RP2_arnoux_yoccoz()
            -a -a -b -b -c -c
            Moebius band
            Lengths: (0.209821688804, 0.114077746827, 0.176100564369)

        """
        sf = arnoux_yoccoz_factor(3)
        return Foliation('a a b b c c', 'moebius',
                [1/sf + 1/sf**2, 1/sf**2 + 1/sf**3, 1/sf + 1/sf**3],
                flips = 'abc')

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
        
        return self.is_bottom_side_moebius() == other.is_bottom_side_moebius()\
                and self._gen_perm == other._gen_perm and \
                abs(self._length_twist_vector -
                    other._length_twist_vector) < epsilon

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
        from foliation_latex import FoliationLatex
        return FoliationLatex(self).tikz_picture()



                

    def in_which_interval(self, value, side):
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
        interval = bisect(self._divvalues[side], value)
        if side == 0 or not self.is_bottom_side_moebius():
            if interval == 0:
                to_check = {self._divvalues[side][0]}
            elif interval == len(self._divvalues[side]):
                to_check = {self._divvalues[side][interval - 1],
                        self._divvalues[side][0] + 1}
            else:
                to_check = {self._divvalues[side][k]
                        for k in {interval - 1, interval}}

            if any(abs(value - x) < epsilon for x in to_check):
                raise SaddleConnectionError()
        return self._all_intervals[side][(interval - 1) % 
                self.num_intervals(side)]


    def apply_iet(self, point, containing_interval):
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
        new_int = containing_interval.pair()
        diff = point - containing_interval.endpoint(0)
        if not containing_interval.is_flipped():
            return (mod_one(new_int.endpoint(0) + diff), new_int)
        else:
            return (mod_one(new_int.endpoint(1) - diff), new_int)


    def _rotated_gen_perm(self, top_rotation, bottom_rotation):
        """
        Returns an involution where the top and bottom rows
        are rotated cyclically.

        INPUT:

        - ``top`` - an integer, shift the top letters
          cyclically by this amount

        - ``bottom`` - an integer, shift the bottom letters
          cyclically by this amount

        OUTPUT:

        - Involution - the rotated Involution

        EXAMPLES::

            sage: i = Involution('a b c b','c a d d', \
                    flips='bc');i
            a -b -c -b
            -c a d d
            sage: i.rotated(1, 1)
            -b a -b -c
            d -c a d
            sage: i.rotated(-6, 2)
            -c -b a -b
            d d -c a
            
        """
        from collections import deque
        labels = [deque(self.labels()[side]) for side in {0, 1}]
        rotations = (top_rotation, bottom_rotation)
        for side in {0, 1}:
            labels[side].rotate(rotations[side])
        return GeneralizedPermutation(list(labels[0]), list(labels[1]), 
                flips = self.flips())

    def _reversed_gen_perm(self):
        """
        Returns an involution where the top and bottom rows
        are reversed.

        OUTPUT:

        - Involution - the reversed Involution

        EXAMPLES::

            sage: i = Involution('a b c b','c a d d', \
                    flips='bc');i
            a -b -c -b
            -c a d d
            sage: i.reversed()
            -b -c -b a
            d d a -c

        """
        from collections import deque
        labels = [deque(self.labels()[side]) for side in {0, 1}]
        for side in {0, 1}:
            labels[side].reverse()
        return GeneralizedPermutation(list(labels[0]), list(labels[1]), 
                flips = self.flips())








    DistanceData = namedtuple('DistanceData', 'side, distance, '
            'is_flipped, orig_side, orig_pos')

    TransitionData = namedtuple('TransitionData', 'tr_matrix,new_inv')

    def new_foliation(self, transition_data):
        new_vector = transition_data.tr_matrix * self._length_twist_vector 
        if self._involution.is_bottom_side_moebius():
            return Foliation(transition_data.new_inv, 
                new_vector.list())
        else:
            return Foliation(transition_data.new_inv, 
                new_vector[:-1].list(), new_vector[-1])

    def _simple_transformation_data(self, new_involution, 
            twist_row = None):
        m = matrix(self._num_coeffs)
        for letter in self._involution.alphabet():
            m[new_involution.index(letter),
                    self._involution.index(letter)] = 1
        if not self._involution.is_bottom_side_moebius():
            m[-1] = twist_row
        return self.TransitionData(tr_matrix = m, 
                new_inv = new_involution)


    def _rotation_data(self, k):
        n = len(self._divpoints[0])
        k = k % n
        twist_row = [0] * self._num_coeffs
        twist_row[-1] = 1
        for i in range(k):
            twist_row[self._involution.index(\
                    self._involution[0][n - 1 - i])] += 1
        return self._simple_transformation_data(\
                self._involution.rotated(k, 0),
                twist_row = twist_row)

    def _reverse_data(self):
        return self._simple_transformation_data(\
                self._involution.reversed(),
                twist_row = [0] * (self._num_coeffs - 1) + [-1])

    def rotated(self, k):
        return self.new_foliation(self._rotation_data(k))

    def reversed(self):
        return self.new_foliation(self._reverse_data())


       
    def _create_transition_data(self, intervals, total, 
            distance_getter):
        intersections = self._get_intersections(intervals)

        distances = [distance_getter(intersections, 
            side, pos) for side in {0, 1}
            for pos in range(len(intersections[side]))]

        sorted_distances = self._get_sorted_distances(distances)
        new_involution = self._get_involution(sorted_distances,
                is_lift = False)
        #print intersections
        #print distances
        #print sorted_distances

        m = matrix(self._num_coeffs)
        done = set()
        for side in {0,1}:
            for i in range(len(sorted_distances[side])):
                letter = new_involution[side][i]
                if letter in done:
                    continue
                diff = sorted_distances[side][(i + 1) % 
                        len(sorted_distances[side])].distance -\
                                sorted_distances[side][i].distance
                if diff.value < 0:
                    diff += total
                m[len(done)] = diff.coefficients
                done.add(letter)
        if len(sorted_distances[1]) != 0:
            m[-1] = sorted_distances[1][0].distance.coefficients -\
                    sorted_distances[0][0].distance.coefficients

        return self.TransitionData(tr_matrix = m,
                new_inv = new_involution)

    @classmethod
    def from_separatrices(cls, separatrices, arc_length = 1, is_lift = False):
        foliation = separatrices[0][0].foliation
        done = set()
        flips = set()
        remaining_labels = range((len(separatrices[0]) +
                                  len(separatrices[1]))/2, 0, -1)
        gen_perm = [[0] * len(separatrices[i]) for i in range(2)]
        lengths = {}
        for side in range(2):
            for i in range(len(separatrices[side])):
                if (side, i) in done:
                    continue
                label = remaining_labels.pop()
                end = 0 # traces which end of the current interval we are
                s = separatrices[side][i]
                interval = s.first_interval()
                if s.is_flipped():
                    interval = interval.prev()
                    end = 1
                interval = interval.pair()
                if interval.is_flipped():
                    end = (end + 1) % 2
                if end == 1:
                    interval = interval.next()
                (new_side, new_i) = next((side, j) for side in range(2)
                                         for j in range(len(separatrices[side]))
                                         if interval == separatrices[side][j].first_interval() and
                                         (not is_lift or separatrices[side][j].is_flipped()
                                          == (end==1)))
                
                if separatrices[new_side][new_i].is_flipped():
                    end = (end + 1) % 2
                if end == 1:
                    new_i = (new_i - 1) % len(separatrices[new_side])
                    flips.add(label)
                gen_perm[side][i] = gen_perm[new_side][new_i] = label
                done.add((side, i))
                done.add((new_side, new_i))
            
                s1 = separatrices[side][i]
                s2 = separatrices[side][(i+1)%len(separatrices[side])]
                lengths[label] = mod_one(s2.endpoint() - s1.endpoint())
                if s1.end_side() < s2.end_side():
                    lengths[label] -= 1 - arc_length
    
        if gen_perm[1] == []:
            gen_perm[1] = 'moebius'
            twist = None
        else:
            twist = mod_one(separatrices[1][0].endpoint() -
                            separatrices[0][0].endpoint())
        return Foliation(*gen_perm, lengths, flips = flips, twist = twist)



    def foliation_orientable_double_cover(self):
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

        """
        if self.is_foliation_orientable():
            return self
        separatrices = []
        
        for interval in self.intervals():
            if interval.is_flipped() and interval.prev() == interval:
                continue
            for i in range(2):
                separatrices.append(Separatrix(self, interval,
                                               number_of_flips_to_stop=i))
            
        return Foliation.from_separatrices(Separatrix.sorted_separatrices(
            separatrices), is_lift = True)
 

    def surface_orientable_double_cover(self):
        """
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
        if self._involution.is_surface_orientable():
            return self
        assert(self._involution.is_bottom_side_moebius())
        n = len(self._divpoints[0])
        alphabet = range(n, 0, -1)
        inv_list = [[0] * n for i in range(2)]
        done = set()
        flips = set()
        lengths = {}
        for i in range(n):
            if i in done:
                continue
            j = self._involution.pair((0, i))[1]
            letter = alphabet.pop()
            letter2 = alphabet.pop()
            lengths[letter] = lengths[letter2] = \
                    self._lengths[self._involution[0][i]].value
            if self._involution.is_flipped((0, i)):
                inv_list[0][i] = inv_list[0][j] = letter
                inv_list[1][i] = inv_list[1][j] = letter2
                flips.add(letter)
                flips.add(letter2)
            else:
                inv_list[0][i] = inv_list[1][j] = letter
                inv_list[0][j] = inv_list[1][i] = letter2
            done.add(i); done.add(j)
        return Foliation(Involution(*inv_list, flips = flips), 
                lengths, twist = Rational('1/2'))



    def _check_not_flipped(self, side, pos):
        if self._involution.is_flipped(side, pos):
            raise RestrictionError('The specified interval should '
                    'not be flipped, but it is.')

    def _check_flipped(self, side, pos):
        if not self._involution.is_flipped(side, pos):
            raise RestrictionError('The specified interval '
                    'should be flipped, but it isn\'t.')

    def _check_same_side(self, side, pos):
        if self._involution.pair((side, pos))[0] != side:
            raise RestrictionError('The interval pair should '
                    'be on the same same side, but it isn\'t.')

    def _check_different_side(self, side, pos):
        if self._involution.pair((side, pos))[0] == side:
            raise RestrictionError('The interval pair should '
                    'be on different sides, but it isn\'t.')

    def _assert_different_pair(self, side1, pos1, side2, pos2):
        if (side1, pos1) in {(side2, pos2),
                self._involution.pair((side2, pos2))}:
            raise RestrictionError('The chosen pairs should be '
                    'different, but they are not.')

    def _get_intersections(self, intervals = None):
        return [[self._first_intersection(side, pos, 
            intervals = intervals)
                for pos in range(len(self._divpoints[side]))]
                for side in {0,1}]

    def _get_right_endpoint(self, side, pos):
        return self._divpoints[side][(pos + 1) %
                len(self._divpoints[side])]

    def _get_left_endpoint(self, side, pos):
        return self._divpoints[side][pos]


    def _restrict_not_flipped_same_side(self, side, pos):
        self._check_not_flipped(side, pos)
        self._check_same_side(side, pos)
        #n = len(self._divpoints[side])
        left_endpoint = self._get_right_endpoint(side, pos)
        side2, pos2 = self._involution.pair((side, pos))
        right_endpoint = self._get_right_endpoint(side2, pos2)

        interval = Arc(left_endpoint, right_endpoint)
        total = interval.length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            distance = (int_data.point - left_endpoint).mod_one()
            if int_data.side == 1:
                distance += total
            return self.DistanceData(side = 0,
                distance = distance,
                is_flipped = int_data.is_flipped,
                orig_side = side,
                orig_pos = pos) 

        return self._create_transition_data([interval], total + total,
                get_distance_data)


    def _restrict_not_flipped_different_side(self, side, pos):
        self._check_not_flipped(side, pos)
        self._check_different_side(side, pos)
        
        left_endpoint = self._get_right_endpoint(side, pos)
        side2, pos2 = self._involution.pair((side, pair))
        right_endpoint = self._get_right_endpoint(side2, pos2)

        interval = Arc(left_endpoint, right_endpoint)

        total = interval.length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            return self.DistanceData(side = int_data.side,
                distance = (int_data.point -left_endpoint).mod_one(),
                is_flipped = int_data.is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data([interval], total)

           
    def _restrict_flipped_two_sided(self, side1, pos1,
            side2, pos2):
        self._assert_different_pair(side1, pos1, side2, pos2)

        Endpoint = namedtuple('Endpoint', 
                'point, end, is_closed, side')

        ints =[[(side1, pos1),
            self._involution.pair((side1, pos1))], [(side2, pos2), 
                self._involution.pair((side2, pos2))]]
        endpoints = []
        for i in {0,1}:
            self._check_flipped(*ints[i][0])
            endpoints.extend([Endpoint(point=self._get_left_endpoint(\
                    *ints[i][0]),end = 'left', is_closed = (i == 1), 
                    side = ints[i][0][0]), 
                    Endpoint(\
                    point = self._get_right_endpoint(*ints[i][1]),
                    end = 'right', is_closed = (i == 1),
                    side = ints[i][1][0])])

        endpoints.sort(key = lambda x: x.point.value)
        while endpoints[0].end == 'right' or\
                endpoints[1].end == 'left':
                    endpoints = endpoints[-1:] + endpoints[:-1]


        def create_interval(i, j):
            return Arc(endpoints[i].point,
                endpoints[j].point,
                left_openness = endpoints[i].is_closed,
                right_openness = endpoints[j].is_closed)


        if endpoints[1].end == 'right':
            intervals = [create_interval(0, 1), create_interval(2, 3)]
        else:
            if endpoints[1].side == endpoints[2].side:
                interval1 = create_interval(1, 2)
                interval2 = create_interval(0, 3)
            else:
                interval2 = create_interval(0, 2)
                interval1 = create_interval(1, 3)
            if endpoints[1].side == 0:
                intervals = [interval1, interval2]
            else:
                intervals = [interval2, interval1]

        total = intervals[0].length() + intervals[1].length() 

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            if int_data.side == 0:
                if intervals[0].contains(int_data.point):
                    distance = (int_data.point -
                            points_with_coeffs[intervals[0][0]]).\
                                    mod_one()
                    new_side = 0
                    is_flipped = int_data.is_flipped
                else:
                    distance = total - (int_data.point - 
                            points_with_coeffs[intervals[1][0]]).\
                                    mod_one()
                    new_side = 1
                    is_flipped = not int_data.is_flipped
            else:
                if intervals[1].contains(int_data.point):
                    distance = total - (int_data.point - 
                            points_with_coeffs[intervals[1][0]]).\
                                    mod_one()
                    new_side = 0
                    is_flipped = not int_data.is_flipped
                else:
                    distance = (int_data.point -
                            points_with_coeffs[intervals[0][0]]).\
                                    mod_one()
                    new_side = 1
                    is_flipped = int_data.is_flipped
            return self.DistanceData(side = new_side,
                distance = distance,
                is_flipped = is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data(intervals, total,
                get_distance_data)

    @staticmethod
    def _get_distances(intersections, getter_function):
        return [getter_function(side, pos) for side in {0, 1}
            for pos in range(len(intersections[side]))]



    def _restrict_flipped_same_side_one_sided(self, side1, pos1,
            side2, pos2):
        self._check_flipped(side1, pos1)
        self._check_flipped(side2, pos2)
        self._assert_different_pair(side1, pos1, side2, pos2)
        assert(self._involution.is_bottom_side_moebius())
        
        n = len(self._divpoints[0])
        center_left = self._divpoints[0][(pos1 + 1) % n]
        center_right = self._divpoints[0][pos2]
        other_left = self._divpoints[0][self._involution.pair(\
                (side1, pos1))[1]]
        other_right = self._divpoints[0][(self._involution.pair(\
                (side2, pos2))[1] + 1) % n]
        center_int = Arc(center_left, center_right,
                left_openness = True, right_openness = True)
        if center_int.contains(other_left) or \
                center_int.contains(other_right):
                    raise RestrictionError('Specified transverse '
                    'curve does not exist. Combinatorically invalid '
                    'choice for the intervals.')
        if (other_right - other_left).mod_one().value <= Rational('1/2'):
            raise RestrictionError('Specified transverse curve does '
            'not exist. The length parameters doesn\'t allow '
            'closing the curve to a transverse curve.')
        other_left = other_left.half_added()

        other_int = Arc(other_left, other_right,
                left_openness = False, right_openness = False)

        intervals = [center_int, other_int]

        total = intervals[0].length() + intervals[1].length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            if int_data.side == 0:
                if intervals[0].contains(int_data.point):
                    distance = (int_data.point - center_left).\
                            mod_one()
                    is_flipped = int_data.is_flipped
                else:
                    distance = total + total - (int_data.point - 
                            other_left).mod_one()
                    is_flipped = not int_data.is_flipped
            else:
                if intervals[1].contains(int_data.point):
                    distance = total - (int_data.point - 
                            other_left).mod_one()
                    is_flipped = not int_data.is_flipped
                else:
                    distance = (int_data.point - center_left).\
                            mod_one() + total
                    is_flipped = int_data.is_flipped
            return self.DistanceData(side = 0,
                distance = distance,
                is_flipped = is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data(intervals, total + total,
                get_distance_data)



#    def _find_pseudo_anosov_candidates(self, depth,


#f = Foliation(Involution('a a b b c c'), [1, 2, 3])

