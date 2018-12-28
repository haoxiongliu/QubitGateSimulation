import numpy as np


def gate_gene(name, n, k, c=-1):
    cell_H = np.array([[2 ** -0.5, 2 ** -0.5], [2 ** -0.5, -(2 ** -0.5)]], dtype=complex)
    cell_X = np.array([[0, 1], [1, 0]], dtype=complex)
    cell_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    cell_Z = np.array([[1, 0], [0, -1]], dtype=complex)
    cell_T = np.array([[1, 0], [0, (2 ** -0.5 + (2 ** -0.5) * 1j)]], dtype=complex)
    cell_S = cell_T.dot(cell_T)
    cell_T_her = np.conjugate(np.transpose(cell_T))
    dic = dict(H=cell_H, T=cell_T, X=cell_X, Y=cell_Y, Z=cell_Z, S=cell_S,
               T_her=cell_T_her)
    cell = dic[name]

    N = 2 ** n
    a = np.identity(N, dtype=complex)
    if c == k:
        raise ValueError
    if c in range(n):
        con_range = [m * 2 ** (n - c) + 2 ** (n - c - 1) + i for m in range(2 ** c) for i in range(2 ** (n - c - 1))]
    elif c == -1:
        con_range = range(2 ** n)
    else:
        raise ValueError
    for i_0 in con_range:
        if (i_0 // (2 ** (n - k - 1))) % 2 == 0:
            i_1 = i_0 + 2 ** (n - k - 1)
            a[i_0, i_0] = cell[0, 0]
            a[i_0, i_1] = cell[0, 1]
            a[i_1, i_0] = cell[1, 0]
            a[i_1, i_1] = cell[1, 1]
    return a
