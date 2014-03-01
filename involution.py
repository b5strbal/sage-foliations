from sage.structure.sage_object import SageObject

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








    



