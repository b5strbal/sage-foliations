from sage.structure.sage_object import SageObject
from collections import deque

class OrientedGraph(SageObject):
    def __init__(self, square_matrix):
        self._n = square_matrix.nrows()
        self._adj_list = matrix_to_adj_list(square_matrix)
        self._backward_adj_list = matrix_to_adj_list(square_matrix.transpose())
            

    def bfs_distances(self, v, backward = False):
        adj_list = self._backward_adj_list if backward else self._adj_list
        done = set([v])
        queue = deque([v])
        distances = [None] * self._n
        distances[v] = 0
        while len(queue) > 0:
            u = queue.popleft()
            for w in adj_list[u]:
                if w not in done:
                    queue.append(w)
                    done.add(w)
                    distances[w] = distances[u] + 1
        return distances

    def is_primitive(self):
        P = self.bfs_distances(0, backward = True)
        Q = self.bfs_distances(0, backward = False)
        if None in P or None in Q: # not strongly connected
            return False
        C = {Q[i] + 1 + P[j] for i in xrange(self._n)
             for j in self._adj_list[i]} # cycle lengths including 0
        # checking for self-loops
        if any(i in self._adj_list[i] for i in xrange(self._n)):
            return True
        # once it is strongly connected, it is primitive if and only if
        # the gcd of cycle lengths in 1
        return gcd(list(C)) == 1
            
def matrix_to_adj_list(square_matrix):
    n = square_matrix.nrows()
    return [[j for j in xrange(n) if square_matrix[i,j] > 0]
            for i in xrange(n)]
        


def is_perron_frobenius(square_matrix):
    g = OrientedGraph(square_matrix)
    return g.is_primitive()


def pf_eigen_data(square_matrix, field = RDF):
    m = matrix(square_matrix, field)
    evr = m.eigenvectors_right()
    largest = max(evr, key = lambda x: abs(x[0]))
    v = largest[1][0]
    v /= sum(v)
    return (largest[0], v)
    






# EPS = 1e-15
# def pf_eigen_data(square_matrix):
#     m = square_matrix
#     # this n^2 - 2n + 2 is large enough power to make all entries positive,
#     # although it is not important here
#     M = power_up(m)
    
#     # get estimate for pf_eigenvalue and normalize the matrix
#     # ev = sum(M*M.column(0)) / float(sum(M.column(0)))
#     # print ev
#     # print M/ev
#     # print (M/ev).eigenvalues()

#     v = M.column(0)
#     v /= float(sum(v))
#     count = 0
#     while True:
#         count += 1
#         old_v = v
#         # print v
#         v = M*v
#         v /= sum(v)
#         if abs(old_v - v) < EPS:
#             break
#         if count > 100:
#             raise RuntimeError("Eigenvector doesn't seem to converge")
#     return (sum(m*v)/sum(v), v)
            

# LIMIT = 10**30
# def power_up(m):
#     while sum(m.list()) < LIMIT:
#         m = m*m
#     return m
