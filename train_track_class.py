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



    





    
