import matplotlib.pyplot as plt
from pylab import mpl

from clsQubit import *
from gatelib import *

mpl.rcParams['font.sans-serif'] = ['FangSong']  # 指定默认字体
mpl.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题

n = 4

# Construct the matrix U
phi = np.array([int(0b010) / (2 ** 3), int(0b000111) / (2 ** 5)], dtype=complex)
chara = np.exp(2 * np.pi * phi * 1j)
chara_vec = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
mat_U = np.zeros((2, 2), dtype=complex)
for i in range(2):
    mat_U += chara[i] * np.kron(chara_vec[i], chara_vec[i]).reshape((2, 2))

init_super = product_all((('H', n, 0), ('H', n, 1), ('H', n, 2)))
mat_phase = product_all((('typein', n, 3, 2, None, 0, power_mat(mat_U, 1)),
                         ('typein', n, 3, 1, None, 0, power_mat(mat_U, 2)),
                         ('typein', n, 3, 0, None, 0, power_mat(mat_U, 4))))
out_gate = product_all((('H', n, 0), ('R_her', n, 1, 0, None, 2), ('H', n, 1),
                        ('R_her', n, 2, 0, None, 3), ('R_her', n, 2, 1, None, 2),
                        ('H', n, 2), ('H', n, 3)))

# Construct the circuit
cir_U = init_super
cir_U = mat_phase.dot(cir_U)
cir_U = out_gate.dot(cir_U)

# Initial qubits
q_vec = QubitVec(n)
init = np.zeros(2 ** n)
init[0] = 1 / np.sqrt(2)
init[1] = 1 / np.sqrt(2)
q_vec.init_state('typein', init)

q_vec.in_gate(cir_U)
q_vec.calc_pos()
print(cir_U)
print(q_vec.pos_1)
print(q_vec.calc_pos_state((0, 1, 0, 0)))

x = np.linspace(0, 7, 100)
y = x
plt.plot(x, y, linestyle='--')
plt.plot([i for i in range(1, 7)], [i for i in range(1, 7)], 'ro', markersize=5, label=u'模拟数据')
plt.ylabel(u"估计精确比特数$m$")
plt.xlabel(u"量子比特数$n$")
plt.title(u"相位估计算法精确比特位数$m$与量子比特数$n$关系")
plt.text(4, 3, '$m = n$')
plt.axis([0, 7, 0, 7])
plt.show()

#
# p_theo = (np.cos(np.pi*phi[0]))**2
# accur = abs(np.log2(abs(p_estimate - p_theo)))
#
#
# x = np.linspace(1, N_max, 1000)
# y = 0.75*np.log2(x)
# plt.plot(x, y, linestyle='--')
# plt.plot([i for i in range(0, N_max)], accur, 'ro', markersize=1, label=u'模拟数据')
# plt.ylabel(u"估计精确比特数m")
# plt.xlabel(u"迭代次数N")
# plt.title(u"Kitaev算法精确比特位数m与迭代次数N关系")
# plt.text(150, 6, '$m = clog_{2}N$')
# plt.axis([0, 200, 0, 8])
# plt.legend()
# plt.show()
