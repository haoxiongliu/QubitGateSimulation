import math

import numpy as np


def gate_gene(name, n, k, con_q=None, neg_q=None, n_of_R=0, cell_ele=None):
    cell_H = np.array([[2 ** -0.5, 2 ** -0.5], [2 ** -0.5, -(2 ** -0.5)]], dtype=complex)
    cell_X = np.array([[0, 1], [1, 0]], dtype=complex)
    cell_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    cell_Z = np.array([[1, 0], [0, -1]], dtype=complex)
    cell_T = np.array([[1, 0], [0, (2 ** -0.5 + (2 ** -0.5) * 1j)]], dtype=complex)
    cell_S = cell_T.dot(cell_T)
    cell_T_her = np.conjugate(np.transpose(cell_T))
    cell_R = np.array([[1, 0], [0, np.exp(2 * math.pi * 1j / 2 ** n_of_R)]], dtype=complex)
    cell_R_her = np.conjugate(np.transpose(cell_R))
    dic = dict(H=cell_H, T=cell_T, X=cell_X, Y=cell_Y, Z=cell_Z, S=cell_S,
               T_her=cell_T_her, R=cell_R, R_her=cell_R_her)

    if name == 'typein':
        if type(k) == int:
            cell = np.array(cell_ele, dtype=complex).reshape((2, 2))
        else:
            cell = np.array(cell_ele, dtype=complex).reshape((2 ** len(k), 2 ** len(k)))
    else:
        cell = dic[name]

    N = 2 ** n
    a = np.identity(N, dtype=complex)
    con_range = range(N)
    if type(con_q) == int:
        con_q = (con_q,)

    if type(con_q) == tuple:
        for c in con_q:
            if c in range(n) and c != k:
                con_range_temp = [m * 2 ** (n - c) + 2 ** (n - c - 1) + i
                                  for m in range(2 ** c) for i in range(2 ** (n - c - 1))]
                con_range = list(set(con_range).intersection(set(con_range_temp)))
            else:
                raise ValueError

    if type(neg_q) == int:
        neg_q = (neg_q,)

    if type(neg_q) == tuple:
        for c in neg_q:
            if c in range(n) and c != k:
                con_range_temp = [m * 2 ** (n - c) + 2 ** (n - c - 1) + i
                                  for m in range(2 ** c) for i in range(2 ** (n - c - 1))]
                con_range = list(set(con_range).difference(set(con_range_temp)))
            else:
                raise ValueError

    if type(k) == int:
        for i_0 in con_range:
            if (i_0 // (2 ** (n - k - 1))) % 2 == 0:
                i_1 = i_0 + 2 ** (n - k - 1)
                a[i_0, i_0] = cell[0, 0]
                a[i_0, i_1] = cell[0, 1]
                a[i_1, i_0] = cell[1, 0]
                a[i_1, i_1] = cell[1, 1]
            else:
                pass
        return a
    else:
        for i_0 in con_range:
            flag = 1
            for ele_k in k:
                if (i_0 // (2 ** (n - ele_k - 1))) % 2 == 1:
                    flag = 0
                    break
            if flag == 1:
                # index_k = i_0*np.ones(2**len(k))
                index_k = [i_0, i_0, i_0, i_0]
                for i in range(2 ** len(k)):
                    for l in range(len(k)):
                        index_k[i] += 2 ** (n - k[l] - 1) * ((i // 2 ** (len(k) - l - 1)) % 2)
                for i_1 in range(2 ** len(k)):
                    for i_2 in range(2 ** len(k)):
                        i_1a = index_k[i_1]
                        i_2a = index_k[i_2]
                        a[i_1a, i_2a] = cell[i_1, i_2]
            else:
                pass
        return a


def product_all(g_input):
    result = gate_gene(*g_input[0])
    for i in range(1, len(g_input)):
        result = gate_gene(*g_input[i]).dot(result)
    return result


def power_mat(mat, k):
    result = mat
    if k > 1:
        for i in range(k - 1):
            result = result.dot(mat)
    return result
