from separatrix import Separatrix


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
    def __init__(self, sep1, openness1, sep2, openness2):
        self._arc = Arc(sep1.endpoint, sep2.endpoint, openness1, openness2)
        self._sep = [sep1, sep2]
        self._direction = 'right' if openness1 == 'closed' else 'left'

    def _repr_(self):
        s = ''
        s += 'Arc: ' + repr(self._arc) + '(' + repr(self._arc.length()) + ')'
        return s

    def coding(self):
        return (self._sep[0])

    @classmethod
    def get_transverse_curves(cls, sep1, end1, sep2, end2):
        r"""
        
        INPUT:

        - ``sep1`` -- 

        - ``sep2`` -- 

        """
        sep = [sep1, sep2]
        ends = [end1, end2]
        # Checking if the curve is transverse
        if (sep1.is_flipped() != sep2.is_flipped() != \
                sep1.first_interval(end1).is_flipped()):
            return []
            
        if sep1.is_flipped() == (end1 == 1):
            openness = ('open','closed')
        else:
            openness = ('closed','open')
                           
        # Checking if one of the arcs is simple        
        arcs = [Arc(sep[i].endpoint, sep[(i+1)%2].endpoint,
                    *openness) for i in range(2)]
        good_indices = [0,1]
        for i in xrange(2):
            # if the curve is one-sided and the new length after cutting
            # is more than 1, we can't compute the train track transition map
            if arcs[i].length() > 0.5 and sep1.end_side == sep2.end_side:
                good_indices.remove(i)

        import itertools
        for x in itertools.chain(sep1.intersections()[:-1], sep2.intersections()[:-1]):
            for i in good_indices:
                if arcs[i].contains(x):
                    good_indices.remove(i)
                    if len(good_indices) == 0:
                        return []
                    break
                        
        return [TransverseCurve(sep[i], openness[0], sep[(i+1)%2],
                                openness[1]) for
                i in good_indices]
                           
    def is_one_sided(self):
        return self._sep[0].end_side == self._sep[1].end_side

    def new_foliation(self):
        foliation = self._sep[0].foliation
        separatrices = Separatrix.get_all(foliation, self._arc)
               
        from transition_map import new_foliation
        return new_foliation(separatrices, self._sep[0].endpoint,
                             self._sep[0].start_side,
                             is_one_sided = self.is_one_sided(),
                             direction = self._direction,
                             ending_point = self._sep[1].endpoint)

Coding = namedtuple("Coding", "side, index, end")

def get_transverse_curve(foliation, coding):
    interval = foliation.interval(coding.side, coding.index)
    if not interval.is_flipped():
        sep1 = Separatrix(foliation, interval,
                          number_of_flips_to_stop = 0)
        sep2 = Separatrix(foliation, interval.pair(),
                          number_of_flips_to_stop = 0)
        return TransverseCurve.get_transverse_curves(\
                                            sep1, 0, sep2, 0))
    other_int = interval.pair()
    sep1 = Separatrix(foliation, interval,
                      number_of_flips_to_stop = 0)
    sep2 = Separatrix(foliation, other_int.next(),
                      number_of_flips_to_stop = 1)
    curves.extend(TransverseCurve.get_transverse_curves(\
                                        sep1, 0, sep2, 1))

    sep1 = Separatrix(foliation, interval.next(),
                      number_of_flips_to_stop = 0)
    sep2 = Separatrix(foliation, other_int,
                      number_of_flips_to_stop = 1)
    curves.extend(TransverseCurve.get_transverse_curves(\
                                    sep1, 1, sep2, 0))

def get_transverse_curves(foliation):
    curves = []
    done = set()
    for interval in foliation.intervals():
        if not interval.is_flipped():
            if interval in done:
                continue
            done.add(interval.pair())
            sep1 = Separatrix(foliation, interval,
                              number_of_flips_to_stop = 0)
            sep2 = Separatrix(foliation, interval.pair(),
                              number_of_flips_to_stop = 0)
            curves.extend(TransverseCurve.get_transverse_curves(\
                                                sep1, 0, sep2, 0))
        else:
            other_int = interval.pair()
            sep1 = Separatrix(foliation, interval,
                              number_of_flips_to_stop = 0)
            sep2 = Separatrix(foliation, other_int.next(),
                              number_of_flips_to_stop = 1)
            curves.extend(TransverseCurve.get_transverse_curves(\
                                                sep1, 0, sep2, 1))

            sep1 = Separatrix(foliation, interval.next(),
                              number_of_flips_to_stop = 0)
            sep2 = Separatrix(foliation, other_int,
                              number_of_flips_to_stop = 1)
            curves.extend(TransverseCurve.get_transverse_curves(\
                                            sep1, 1, sep2, 0))
            

    return curves

# PseudoAnosov = namedtuple("PseudoAnosov", "start_fol, 

def find_pseudo_anosovs(foliation, depth, tt_map_so_far, coding_so_far):
    curves = get_transverse_curves(foliation)
    
