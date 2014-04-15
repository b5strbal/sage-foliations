from foliation import Foliation
from transverse_curve import tt_map_from_codings, Coding
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import \
    PolynomialRing

# charpoly: x^6-x^5-kx^3-x-1
f = Foliation('1 2 3 4 5 6','5 2 1 4 3 6',[0.170886909342, 0.189075941817, 0.108600690819, 0.0879378482553, 0.0369647624731, 0.406533847294], twist = 0.0914086480774)


def family_A_factor(n, k, field = RDF):
    r"""
    Returns the Perron root of the polynomials 
    $x^n - kx^{n-1} - \cdots - kx - 1$.

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
    R = PolynomialRing(ZZ, 't')
    poly = R([-1] + [-k] * (n-1) + [1])
    return max([abs(x[0]) for x in poly.roots(field)])
    



def family_A_foliation_RP2(k):
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
    sf = family_A_factor(3, k)
    a, b, c = [1/sf, 1/sf**2, 1 - 1/sf - 1/sf**2]
    return Foliation('a a b b c c', 'moebius',
                     [a + b, b + c, c + a],
                     flips = 'abc')



def family_A_foliation(n, k, is_orientable):
    sf = family_A_factor(n, k)
    lengths = [1/sf**i for i in range(1,n)]
    lengths.append(1-sum(lengths))

    if is_orientable:
        top = range(1, 2*n + 1)
        def switch(k):
            if k % 2 == 0:
                return k + 1
            return k - 1
        bottom = [top[switch(i)] for i in range(2*n)]
        new_list = []
        for x in lengths:
            for i in xrange(2):
                new_list.append(x)
        return Foliation(top, bottom, new_list, twist = 1)
    else:
        top = sorted(2 * range(1, n + 1))
        return Foliation(top, 'moebius', lengths)




def orientableAY_tt_map():
    fol = Foliation.orientable_arnoux_yoccoz(3)
    coding_list = [Coding(1,4,0,0,0),Coding(1,1,0,0,0),False,5]
    return tt_map_from_codings(fol, coding_list)


def family_A_tt_map(n, k, is_orientable):
    fol = family_A_foliation(n, k, is_orientable)
    if not is_orientable:
        coding_list = [Coding(0,0,0,0,0),True,-1]
    else:
        b = 2*n-4 if k > 1 else 1
        coding_list = [Coding(0,0,0,0,0),Coding(1,b,0,0,0),False,-1]

    return tt_map_from_codings(fol, coding_list)

