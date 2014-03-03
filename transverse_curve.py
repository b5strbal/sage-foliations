class Arc(SageObject):
    """
    An interval of the unit interval $[0,1]$ with opposite sides 
    identified.

    INPUT:

    - ``left_endpoint`` - PointWithCoefficients
    - ``right_endpoint`` - PointWithCoefficients
    - ``left_openness`` - 'closed' or 'open'. The default is 'closed'
    - ``right_openness`` - 'closed' or 'open'. The default is 'open'

    EXAMPLES:

    Any usual interval can be specified::

        sage: from sage.dynamics.foliations.foliation import Arc, PointWithCoefficients
        sage: Arc(PointWithCoefficients(1/2, [2, 1]), \
                PointWithCoefficients(3/5, [4, -1]))
        [(1/2, (2, 1)), (3/5, (4, -1)))
        sage: Arc(PointWithCoefficients(0, [0]), \
                PointWithCoefficients(0.4, [1]), 'open', 'closed')
        ((0, (0)), (0.400000000000000, (1))]
        
    But intervals wrapping around the endpoints are also considered::

        sage: Arc(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        [(1/2, (0, 1)), (1/4, (1, -2))]

    One can get the endpoints by indexing and the length() using the 
    length()::

        sage: i = Arc(PointWithCoefficients(1/5, (1, 2)), \
                PointWithCoefficients(4/7, (6, -1)))
        sage: i[0]
        (1/5, (1, 2))
        sage: i[1]
        (4/7, (6, -1))
        sage: i.length()
        (13/35, (5, -3))

        sage: j = Arc(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        sage: j.length()
        (3/4, (2, -3))
    
    """
    def __init__(self, left_endpoint, right_endpoint,
            left_openness = 'closed', right_openness = 'open'):
        if not 0 <= left_endpoint < 1 or not \
                0 <= right_endpoint < 1:
                    raise ValueError("The endpoints of the Arc "
                            "must be between 0 and 1.")
        self._lep = left_endpoint
        self._rep = right_endpoint
        if not {'closed', 'open'} >= {left_openness, right_openness}:
            raise ValueError('Openness arguments should be either '
                    '\'closed\' or \'open\'')
        self._left_openness = left_openness
        self._right_openness = right_openness

    @staticmethod
    def _less(x, y, is_equality_allowed):
        """
        Unifies the < and <= operators in a single function.

        INPUT:

        - ``x`` - a number
        - ``y`` - another number
        - ``is_equality_allowed`` - if True, the function returns 
          x<=y, otherwise x<y

        TESTS::

            sage: from sage.dynamics.foliations.foliation import Arc
            sage: Arc._less(1, 2, True)
            True
            sage: Arc._less(2, 1, True)
            False
            sage: Arc._less(1, 1, True)
            True
            sage: Arc._less(1, 1, False)
            False

        """
        if is_equality_allowed:
            return x <= y
        else:
            return x < y

    def _repr_(self):
        """
        Returns the representation of self.

        TESTS::

        sage: from sage.dynamics.foliations.foliation import Arc, PointWithCoefficients
        sage: Arc(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        [(1/2, (0, 1)), (1/4, (1, -2))]

        """
        s = ''
        if self._left_openness == 'closed':
            s += '['
        else:
            s += '('
        s += str(self._lep) + ', ' + str(self._rep)
        if self._right_openness == 'closed':
            s += ']'
        else:
            s += ')'
        return s

    def __getitem__(self, index):
        """
        Returns the left or right endpoint of the interval.

        INPUT:

        - ``index`` - either 0 or 1, for the left and right endpoints
          respectively

        OUTPUT:

        - PointWithCoefficients - one of the endpoints

        EXAMPLES::

            sage: from sage.dynamics.foliations.foliation import Arc, PointWithCoefficients
            sage: i = Arc(PointWithCoefficients(1/5, (1, 2)), \
                    PointWithCoefficients(4/7, (6, -1)))
            sage: i[0]
            (1/5, (1, 2))
            sage: i[1]
            (4/7, (6, -1))

        """
        if index == 0:
            return self._lep
        if index == 1:
            return self._rep

    def length(self):
        """
        Returns the length of the interval.

        OUTPUT:

        - PointWithCoefficients - the length of the interval

        EXMAPLES::

            sage: from sage.dynamics.foliations.foliation import Arc, PointWithCoefficients
            sage: i = Arc(PointWithCoefficients(1/5, (1, 2)), \
                    PointWithCoefficients(4/7, (6, -1)))
            sage: i.length()
            (13/35, (5, -3))

            sage: j = Arc(PointWithCoefficients(0.5, [0, 1]), \
                PointWithCoefficients(0.25,[1, -2]),'closed','closed')
            sage: j.length()
            (0.750000000000000, (2, -3))
    
        """
        return mod_one(self[1] - self[0])
        #return (self._rep - self._lep).mod_one()

    def contains(self, point):
        """
        Decides if a point in contained in self.

        The cooefficients of the point don't matter, only the value.

        INPUT:

        - ``point`` - PointWithCoefficients or a real number, 
            must be in the unit interval $[0, 1)$ 

        OUTPUT:

        - boolean - True if ``point`` is contained in self,
          False if not

        EXAMPLES::

            sage: from sage.dynamics.foliations.foliation import Arc, PointWithCoefficients
            sage: i = Arc(PointWithCoefficients(0.25, [0, 1]), \
                PointWithCoefficients(0.5,[1, -2]),'open','closed')
            sage: i.contains(PointWithCoefficients(0.1, [1,1]))
            False
            sage: i.contains(PointWithCoefficients(0.3, [-2,6]))
            True
            sage: i.contains(PointWithCoefficients(0.25, [-1,5]))
            False
            sage: i.contains(0.25)
            False
            sage: i.contains(PointWithCoefficients(0.5, [1,7]))
            True
            sage: i.contains(0.5)
            True

            sage: j = Arc(PointWithCoefficients(0.5, [0, 1]), \
                PointWithCoefficients(0.25,[1, -2]),'closed','open')
            sage: j.contains(PointWithCoefficients(0.1, [1,1]))
            True
            sage: j.contains(PointWithCoefficients(0.3, [-2,6]))
            False
            sage: j.contains(PointWithCoefficients(0.25, [-1,5]))
            False
            sage: j.contains(PointWithCoefficients(0.5, [1,7]))
            True

            sage: k = Arc(PointWithCoefficients(0.3, (1,1)),\
                    PointWithCoefficients(0.3, (1,1)))
            sage: k.contains(PointWithCoefficients(0.3, (1,1)))
            False

            sage: l = Arc(PointWithCoefficients(0.3, (1,1)),\
                    PointWithCoefficients(0.3, (1,1)), 'closed', \
                    'closed')
            sage: l.contains(PointWithCoefficients(0.3, (1,1)))
            True

        """
        if not 0 <= point < 1:
            raise ValueError("Only points in the unit interval can be"
                    " tested for containment")
        if self[0] <= self[1]:
            if self._less(self[0], point, 
                    is_equality_allowed = (self._left_openness == 
                        'closed')) and\
                    self._less(point, self[1], 
                            is_equality_allowed = 
                            (self._right_openness == 'closed')):
                return True
            else:
                return False
        if self._less(self[1], point, 
                is_equality_allowed = 
                (self._right_openness == 'open')) and\
                self._less(point, self[0], 
                    is_equality_allowed =
                        (self._left_openness == 'open')):
            return False
        else:
            return True



class TransverseCurve(SageObject):
    def __init__(self, sep1, sep2):
        if (sep1.first_interval_end() == 0) == sep1.is_flipped():
            openness = ('closed','open')
        else:
            openness = ('open','closed')
        self._arc = Arc(sep1.endpoint(), sep2.endpoint(), *openness)
        self._sep = [sep1, sep2]


    @classmethod
    def get_transverse_curves(cls, sep1, sep2):
        r"""
        
        INPUT:

        - ``sep1`` -- 

        - ``sep2`` -- 

        """
        sep = [sep1, sep2]
        # Checking if the curve is transverse
        if not (sep1.is_flipped() == sep2.is_flipped() == \
           sep1.traversed_intervals[0].is_flipped()):
            return []
            
        # If the separatrices are short, both arcs are good
        if len(sep1.intersections) == 1 and len(sep2.intersections) == 1:
            return [TransverseCurve(sep1, sep2), TransverseCurve(sep2, sep1)]
                           
        # Checking if one of the arcs is simple        
        arcs = [Arc(sep[i].endpoint(), sep[(i+1)%2].endpoint(),
                    'open', 'open') for i in range(2)] 
        good_indices = [0,1]
        import itertools
        for x in itertools.chain(sep1.intersections, sep2.intersections):
            for i in good_indices:
                if arcs[i].contains(x):
                    good_indices.remove(i)
                    if len(good_indices) == 0:
                        return []
                    break
                        
        return [TransverseCurve(sep[i],sep[(i+1)%2]) for
                i in good_indices]
                           
    def is_one_sided(self):
        return sep1.end_side() == sep2.end_side()

    def cut(self):
        foliation = self._sep[0].foliation
        separatrices = Separatrix.sorted_separatrices(
            Separatrix.get_all(foliation, self._arc), self._arc[0])
            
        if self.is_one_sided():
            separatrices[0].extend(separatrices[1])
            separatrices[1] = []
                
        return new_foliation(separatrices, self._arc.length())

    @classmethod
    def find_transverse_curves(cls, interval, end,
            num_flipped_ints1, num_flipped_ints2, max_lengthening):


                                     
