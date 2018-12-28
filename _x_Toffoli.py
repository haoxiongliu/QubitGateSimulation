import numpy as np

from _x_def_gate import *

# Calculate the final gate
n = 3
# put in gate Toffoli
g_input = (('H', n, 2), ('X', n, 2, 1), ('T_her', n, 2), ('X', n, 2, 0),
           ('T', n, 2), ('X', n, 2, 1), ('T_her', n, 2), ('X', n, 2, 0),
           ('T_her', n, 1), ('T', n, 2), ('X', n, 1, 0), ('H', n, 2),
           ('T_her', n, 1), ('X', n, 1, 0), ('T', n, 0), ('S', n, 1))
cir_U = np.identity(2 ** 3, dtype=complex)
for ele in g_input:
    cir_U = cir_U.dot(gate_gene(*ele))
print(abs(cir_U))

qubit_list = np.array([[q1, q2, q3] for q1 in range(2) for q2 in range(2) for q3 in range(2)])
