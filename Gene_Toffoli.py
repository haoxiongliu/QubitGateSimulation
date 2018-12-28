# generate a toffoli gate with basic gates

import time

from gatelib import *

if __name__ == '__main__':
    n = 3
    # Calculate the final gate
    # put in gates
    g_input = (('H', n, 2), ('X', n, 2, [1]), ('T_her', n, 2), ('X', n, 2, 0),
               ('T', n, 2), ('X', n, 2, 1), ('T_her', n, 2), ('X', n, 2, 0),
               ('T_her', n, 1), ('T', n, 2), ('X', n, 1, 0), ('H', n, 2),
               ('T_her', n, 1), ('X', n, 1, 0), ('T', n, 0), ('S', n, 1))

    start = time.time()
    cir_U = product_all(g_input)
    # cir_U = np.identity(2**n, dtype=complex)
    # for ele in g_input:
    #     cir_U = cir_U.dot(gate_gene(*ele))
    end = time.time()
    exec_time = end - start
    print(exec_time)
    print(abs(cir_U))
