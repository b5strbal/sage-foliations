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

from sage.structure.sage_object import SageObject



class TrainTrackClass(SageObject):
    def __init__(self, tt):
        self._tt_repr = tt

    def _repr_(self):
        s = "Train Track equivalence class with invariants "
        s += repr(self.invariants())
        return s

    def __hash__(self):
        return hash(self.invariants())

    def __eq__(self, other):
        return isinstance(other, TrainTrackClass) and \
            self.invariants() == other.invariants() and \
            len(self._tt_repr.get_symmetries_from(other._tt_repr)) > 0

    def invariants(self):
        return self._tt_repr.invariants()

    def tt_repr(self):
        return self._tt_repr



    





    
