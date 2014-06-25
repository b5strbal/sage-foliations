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

from foliation import Foliation
from transverse_curve import Coding
from train_track_map import tt_map_from_codings
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF

from sage.rings.polynomial.polynomial_ring_constructor import \
    PolynomialRing
from pseudo_anosov import PseudoAnosov
from constants import *
# charpoly: x^6-x^5-kx^3-x-1
# f = Foliation('1 2 3 4 5 6','5 2 1 4 3 6',[0.170886909342, 0.189075941817, 0.108600690819, 0.0879378482553, 0.0369647624731, 0.406533847294], twist = 0.0914086480774)


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
    
# def family_B_factor(n, k, field = RDF):
#     R = PolynomialRing(ZZ, 't')
#     poly = R([1] + [-k+1] * (n-2) + [-k-1,1])
#     print poly
#     return max([abs(x[0]) for x in poly.roots(field)])
    


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

# def family_B_foliation_nonor(n, k):
#     sf = family_B_factor(n, k)
#     lengths = [1/sf**i for i in range(1,n)]
#     lengths.append(1-sum(lengths))
#     top = sorted(2 * range(1, n + 1))
#     return Foliation(top, 'moebius', lengths)


# def orientableAY_tt_map():
#     fol = Foliation.orientable_arnoux_yoccoz(3)
#     coding_list = [Coding(1,4,0,0,0),Coding(1,1,0,0,0),False,5]
#     return tt_map_from_codings(fol, coding_list)


from transverse_curve import RestrictionCoding, Coding
from foliation import SymmetryCoding
from interval import Interval

def family_A_PA(n, k, is_orientable):
    fol = family_A_foliation(n, k, is_orientable)
    if not is_orientable:
        coding_list = [RestrictionCoding(fol, Coding(0,0,1,0,0)),
                       SymmetryCoding(Interval(0,1), RIGHT)]
    else:
        b = 2*n-4 if k > 1 else 1
        coding_list = [RestrictionCoding(fol,Coding(0,0,0,0,0)),Coding(1,b,0,0,0),SymmetryCoding(Interval(0,1), RIGHT)]

    tt_map = tt_map_from_codings(fol.train_track, coding_list)
    # return (tt_map.domain, tt_map.codomain)
    return PseudoAnosov(tt_map)


fol_min4 = Foliation('e e a b f f d b d c a c', 'moebius', (1, 0.5834783686864043, 0.8827306520664273, 0.6609925318901200, 1.406064340122059, 0.9293980281776905))
tt_map = tt_map_from_codings(fol_min4.train_track,
                             [RestrictionCoding(fol_min4, Coding(0,3,1,0,0))])
pa_min4 = PseudoAnosov(tt_map)


# def family_B_tt_map(n, k):
#     fol = family_B_foliation_nonor(n, k)
#     coding_list = [Coding(0,0,1,0,0),False,1]
#     return tt_map_from_codings(fol, coding_list)


def family_A_PA_other(n, k, is_orientable):
    fol = family_A_foliation(n, k, is_orientable)
    if not is_orientable:
        coding_list = [RestrictionCoding(fol, Coding(0,0,0,0,0)),
                       SymmetryCoding(Interval(0,5), LEFT)]
    else:
        b = 2*n-4 if k > 1 else 1
        coding_list = [RestrictionCoding(fol,Coding(0,0,0,0,0)),Coding(1,b,0,0,0),SymmetryCoding(Interval(0,1), RIGHT)]

    tt_map = tt_map_from_codings(fol.train_track, coding_list)
    # return (tt_map.domain, tt_map.codomain)
    return PseudoAnosov(tt_map)
