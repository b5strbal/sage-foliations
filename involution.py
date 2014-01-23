from collections import deque
from sage.structure.sage_object import SageObject
from sage.dynamics.interval_exchanges.constructors import GeneralizedPermutation

class Involution(SageObject):
    """
    A wrapper for different kinds of Generalized Permutations.

    Foliation objects are based on Involutions which are 
    basically the same as GeneralizedPermutations already 
    in sage. Involution objects have a few extra methods.
   
    INPUT:

    - ``top_letters`` - a string where interval names are
      separated by spaces, or a list is interval names
    - ``bottom_letters`` - a string where interval names are
      separated by spaces, or a list is interval names

    - ``flips`` - a list or set of flipped interval names, or
      a string containing flipped letters in case the names 
      are single characters

    EXAMPLES:

    The top and bottom intervals can be specified as strings,
    the letters (or numbers or words) separated by spaces. The
    flips can be a string in containing the letters of flipped
    intervals::

        sage: i = Involution('a b c', 'c b a', flips = 'ab'); i
        -a -b  c
         c -b -a

    If the names of intervals are not single characters, this
    flips notation doesn't work::

        sage: i = Involution('1 2 33', '33 2 1', flips = '33'); i
        Traceback (most recent call last):
        ...
        TypeError: The flip list is not valid

    In this case the flipped intervals must be presented in
    a list or set so that element testing works::

        sage: i = Involution('1 2 33', '33 2 1', flips = ['33']); i
         1  2 -33
        -33  2  1
        
    The top and bottom intervals can also be listed in a list::

        sage: i = Involution(['a','b','c'],['c','b','a'], 'ab'); i
        -a -b  c
         c -b -a

    If the second argument is omitted, the bottom side of the
    curve is considered as a Moebius band without punctures::

        sage: i = Involution('a a b b c c', flips = 'abc'); i
        -a -a -b -b -c -c
        Moebius band
        sage: i.singularity_type()
        (3, 1, 1, 1)

    It is only be omitting the bottom letters that one gets
    a Moebius band without punctures. The following results in
    a once punctured Moebius band on the bottom::

        sage: i = Involution('a a b b c c', 'd d', flips = 'abc'); i
        -a -a -b -b -c -c
         d  d
        sage: i.singularity_type()
        (3, 2, 1, 1, 1)

    The mapping class group of the closed torus and the 
    once-punctured torus are the same, and since we are interested in
    constructed pseudo-anosovs, we won't complicate the code with
    treating the case of the closed torus separately. But we can 
    represent once-punctured tori as follows::

        sage: i = Involution('a', 'a'); i
        a
        a
        sage: i.singularity_type()
        (2,)

    And here is a twice-punctured torus::

        sage: i = Involution('a b', 'a b'); i
        a b 
        a b
        sage: i.singularity_type()
        (2, 2)

    """
    def __init__(self, top_letters, bottom_letters = None, 
            flips = []):
        if bottom_letters == None: # bottom side is Moebius
            bottom_letters = 'JOKER JOKER'
        self._gen_perm = GeneralizedPermutation(\
                top_letters, bottom_letters, flips = flips)

        # initializing self._pair
        self._pair = {}
        for i in range(2):
            for j in range(len(self[i])):
                for a in range(2):
                    for b in range(len(self[a])):
                        if (i != a or j != b) and \
                                self[i][j] == self[a][b]:
                                    self._pair[(i,j)] = (a,b)

        # initializing self._index
        self._index = {}
        count = 0
        done = set()
        for (i, j) in sorted(self._pair.keys()):
            letter = self[i][j]
            if letter in done:
                continue
            done.add(letter)
            self._index[letter] = count
            count += 1

        # initializing self._singularity_partition
        done = set()
        if self.is_bottom_side_moebius():
            done = {(1,0), (1,1)}
        partition = []
        for (i, j) in self._pair:
            if (i, j) in done:
                continue
            (a, b) = (i, j)
            partition.append([])
            direction = 'left'
            while True:
                if direction == 'left':
                    (a, b) = self._pair[(a, (b - 1) % 
                        len(self[a]))]
                    if not self.is_flipped((a, b)):
                        b = (b + 1) % len(self[a])
                        direction = 'away'
                else:
                    (a, b) = self._pair[(a,b)]
                    if not self.is_flipped((a, b)):
                        direction = 'left'
                    else:
                        b = (b + 1) % len(self[a])
                partition[-1].append((a, b))
                done.add((a, b))
                if (a, b) == (i, j):
                    break

        self._singularity_partition = partition


    def _repr_(self):
        """
        Returns a representation of self.

        EXAMPLES:

        Usually it is the same as the representation of a
        GeneralizedPermutation::

            sage: i = Involution('a a b b', 'c c', flips = 'ab'); i
            -a -a -b -b
             c  c

        It's different when the bottom side is a Moebius band::

            sage: i = Involution('a a b b', flips ='ab'); i
            -a -a -b -b
            Moebius band

        """
        if self.is_bottom_side_moebius():
            return repr(self._gen_perm).split('\n')[0] + \
                    '\nMoebius band'
        return repr(self._gen_perm)

    def __getitem__(self, index):
        """
        Returns the list of top of bottom letters.

        INPUT:

        - ``index`` - 0 for the top letters, 1 for the bottom

        EXAMPLES::

            sage: i = Involution('a a b b','c c',flips = 'b');i
            a a -b -b
            c c
            sage: i[0]
            ['a', 'a', 'b', 'b']
            sage: i[1]
            ['c', 'c']

        """
        return self._gen_perm.list()[index]

    def __eq__(self, other):
        """
        Decides if two Involutions are the same up to renaming
        letters.

        EXAMPLES::
            
            sage: i = Involution('a a b b', 'c c')
            sage: j = Involution('1 1 2 2', '3 3')
            sage: k = Involution('a a b b', 'c c', flips = 'c')
            sage: l = Involution('a a b b')
            sage: i == j
            True
            sage: i == k
            False
            sage: i == l
            False
        """ 
        if self.is_bottom_side_moebius() !=\
                other.is_bottom_side_moebius():
                    return False
        return self._gen_perm == other._gen_perm

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

    def is_flipped(self, pos):
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
        x = self._gen_perm[pos[0]][pos[1]]
        if isinstance(x, tuple):
            # self._gen_perm is a FlippedLabelledPermutationLI
            return x[1] == -1
        #self._gen_perm is a LabelledPermutationIET 
        return False

    def index(self, letter):
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
        return self._index[letter]

    def pair(self, pos):
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
        return self._pair[pos]

    @classmethod
    def orientable_arnoux_yoccoz(self, genus):
        """
        Returns the Involution of the Arnoux-Yoccoz foliations
        on orientable surfaces.

        INPUT:

        - ``genus`` - the genus of the surface

        OUTPUT:

        - Involution

        EXAMPLES::

            sage: i = Involution.orientable_arnoux_yoccoz(3); i
            1 2 3 4 5 6
            2 1 4 3 6 5
            sage: i = Involution.orientable_arnoux_yoccoz(4); i
            1 2 3 4 5 6 7 8
            2 1 4 3 6 5 8 7

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
        return Involution(top, bottom)

    @classmethod
    def nonorientable_arnoux_yoccoz(self, genus):
        """
        Returns the Involution of the Arnoux-Yoccoz foliations
        on non-orientable surfaces.

        Take care with the the genus here: 
        The orienting double cover of a closed genus 
        $g$ non-orientable surface is the closed genus $g-1$
        orientable surface.
        
        INPUT:

        - ``genus`` - the non-orientable genus of the surface

        OUTPUT:

        - Involution -

        EXAMPLES::

            sage: i = Involution.nonorientable_arnoux_yoccoz(4); i
            1 1 2 2 3 3
            Moebius band
            sage: i = Involution.nonorientable_arnoux_yoccoz(5); i
            1 1 2 2 3 3 4 4
            Moebius band

        """
        if genus < 4:
            raise ValueError('The genus of a non-orientable '
                    'Arnoux-Yoccoz surface is at least 4')
        top = sorted(2 * range(1, genus))
        return Involution(top)

    @classmethod
    def RP2_arnoux_yoccoz(self):
        """
        Returns the Involution of the Arnoux-Yoccoz foliation
        on the projective plane.

        OUTPUT:

        - Involution -

        EXAMPLES::

            sage: i = Involution.RP2_arnoux_yoccoz(); i
            -a -a -b -b -c -c
            Moebius band

        """
        return Involution('a a b b c c', flips = 'abc')

    def _top_deque(self):
        """
        Returns the list of top letters as a deque.

        OUTPUT:

        - deque -- the deque of top letters

        TESTS::

            sage: i = Involution('a a b b','c c', flips='ab')
            sage: i._top_deque()
            deque(['a', 'a', 'b', 'b'])

            sage: i = Involution('a a b b','c c')
            sage: i._top_deque()
            deque(['a', 'a', 'b', 'b'])
        """
        return deque(self._gen_perm.list()[0])

    def _bottom_deque(self):
        """
        Returns the list of bottom letters as a deque.

        OUTPUT:

        - deque -- the deque of bottom letters

        TESTS::

            sage: i = Involution('a a b b','c c', flips='ab')
            sage: i._bottom_deque()
            deque(['c', 'c'])

            sage: i = Involution('a a b b','c c')
            sage: i._bottom_deque()
            deque(['c', 'c'])

        """
        return deque(self._gen_perm.list()[1])

    def rotated(self, top, bottom):
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
        top_list = self._top_deque() 
        bottom_list = self._bottom_deque() 
        top_list.rotate(top)
        bottom_list.rotate(bottom)
        return Involution(list(top_list), list(bottom_list),
                self.flips())


    def reversed(self):
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
        top_list = self._top_deque() 
        bottom_list = self._bottom_deque() 
        top_list.reverse()
        bottom_list.reverse()
        return Involution(list(top_list), list(bottom_list),
                self.flips())

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

    def which_singularity(self, pos):
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
        sp = self._singularity_partition
        for i in range(len(sp)):
            if pos in sp[i]:
                return i
        raise ValueError("Invalid singularity specification.")
         
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
        return self[1] == ['JOKER', 'JOKER']

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
        for x in self._pair:
            a = (x[0] == self._pair[x][0])
            b = self.is_flipped(x)
            if a != b:
                return False
        return True

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




