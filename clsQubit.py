import random

import numpy as np


def gene_sing_qubit():
    a = random.random()
    b = (1 - a ** 2) ** 0.5
    return np.array([a, b], dtype=complex)


def gene_n_qubit(n):
    if n == 1:
        result = gene_sing_qubit()
    else:
        result = np.kron(gene_n_qubit(n - 1), gene_sing_qubit())
    return result


class QubitVec:
    def __init__(self, dim):
        self.dim = dim
        self.vec = []
        self.pos_0 = []
        self.pos_1 = []

    def init_state(self, mode='zero', s=None):
        if mode == 'rand':
            self.vec = gene_n_qubit(self.dim)
        elif mode == 'typein':
            if len(s) != 2 ** self.dim:
                raise ValueError
            self.vec = np.array(s, dtype=complex)
        elif mode == 'zero':
            self.vec = np.zeros(2 ** self.dim, dtype=complex)
            self.vec[0] = 1
        else:
            raise NameError('Invalid mode')

    def in_gate(self, gate):
        if gate.shape == (2 ** self.dim, 2 ** self.dim):
            self.vec = gate.dot(self.vec)
        else:
            raise NameError('Invalid Input')

    def calc_pos(self):
        n = self.dim
        self.pos_0 = []
        self.pos_1 = []
        for k in range(n):
            arr_0 = np.array([self.vec[m * 2 ** (n - k) + i] for m in range(2 ** k) for i in range(2 ** (n - k - 1))])
            self.pos_0.append((abs(arr_0) ** 2).sum())
            self.pos_1.append(1 - self.pos_0[k])
            del arr_0

    def calc_pos_state(self, state=None):
        if len(state) != self.dim:
            raise ValueError
        if state is None:
            state = np.zeros(self.dim)

        index_vec = 0
        for i in range(len(state)):
            index_vec += (2 ** i) * state[len(state) - i - 1]

        return abs(self.vec[index_vec]) ** 2

    def measure(self, k):
        r = random.random()
        if 0 <= r < self.pos_1[k]:
            return 1
        elif self.pos_1[k] <= r < 1:
            return 0
        else:
            raise ValueError
