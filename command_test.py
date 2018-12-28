from gatelib import *

n = 4
grover_diffusion = (('H', n, 0), ('H', n, 1), ('H', n, 2), ('X', n, 0),
                    ('X', n, 1), ('X', n, 2), ('H', n, 2),
                    ('X', n, 2, (0, 1)), ('H', n, 2), ('X', n, 0), ('X', n, 1), ('X', n, 2),
                    ('H', n, 0), ('H', n, 1), ('H', n, 2))

gate_dif = np.identity(2 ** n, dtype=complex)
for ele in grover_diffusion:
    gate_dif = gate_dif.dot(gate_gene(*ele))
print(gate_gene('X', n, 3, (0, 1), (2,)))
