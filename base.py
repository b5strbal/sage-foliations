r"""
Constants, functions, etc. used in multiple files throughout the code.

AUTHORS:

- Balazs Strenner (2014-06-16): initial version


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

LEFT = 0
MID = 0.5
RIGHT = 1
START = 0
END = 1
TOP = 0
BOTTOM = 1

SIGNED = 0
UNSIGNED = 1

TRAIN_TRACK_MODULE = 0
COHOMOLOGY = 1

EPSILON = 1e-10


class MyException(Exception):
    """The base exception class for user-defined exceptions."""
    def __init__(self, value =''):
        self.value = value

    def __str__(self):
        return self.value


class SaddleConnectionError(MyException):
    """The error signaling that a saddle connection is found."""
    pass

class RestrictionError(MyException):
    pass



from sage.functions.other import floor 
def mod_one(x):
    """
    Return ``x`` mod 1.

    INPUT:

    - ``x`` -- a real number

    OUTPUT:

    ``x`` mod 1, a number in `[0,1)`.

    TESTS::

        sage: from sage.dynamics.foliations.base import mod_one
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
    return x - floor(x)

