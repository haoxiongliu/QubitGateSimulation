from matplotlib import pyplot as plt

from clsQubit import *
from gatelib import *

plt.rcParams['font.sans-serif'] = ['SimHei'] # 步骤一（替换sans-serif字体）
plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题

n = 2

# Construct the matrix U
phi = np.array([int(0b00011101) / (2 ** 8), int(0b000111) / (2 ** 5)], dtype=complex)
chara = np.exp(2 * np.pi * phi * 1j)  # 选择特征值
# 选择特征矢量（满足正交）
chara_vec = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
# 生成酉算子U
mat_U = np.zeros((2, 2), dtype=complex)
for i in range(2):
    mat_U += chara[i] * np.kron(chara_vec[i], chara_vec[i]).reshape((2, 2))

init_super = product_all((('H', n, 0),))
mat_phase = gate_gene('typein', n, 1, 0, None, 0, tuple(np.ndarray.flatten(mat_U)))
out_gate = product_all((('H', n, 0),))

q_vec = QubitVec(n)
init = (2 ** -0.5, (2 ** -0.5), 0, 0)
q_vec.init_state('typein', init)

N_max = 200
p_estimate = np.zeros(N_max)
for N in range(1, N_max + 1):
    result = np.zeros(N)
    for i in range(N):
        cir_U = init_super
        cir_U = mat_phase.dot(cir_U)
        cir_U = out_gate.dot(cir_U)

        q_vec.init_state('typein', tuple(init))
        q_vec.in_gate(cir_U)
        q_vec.calc_pos()
        result[i] = 1 - q_vec.measure(0)
    p_estimate[N - 1] = np.sum(result) / len(result)

p_theo = (np.cos(np.pi * phi[0])) ** 2
accur = abs(np.log2(abs(p_estimate - p_theo)))

x = np.linspace(1, N_max, 1000)
y = 0.75 * np.log2(x)
plt.plot(x, y, linestyle='--')
plt.plot([i for i in range(0, N_max)], accur, 'ro', markersize=1, label=u'模拟数据')
plt.ylabel(u"估计精确比特位数$m$")
plt.xlabel(u"调用算子$U$次数$N$")
plt.title(u"Kitaev算法精确比特位数$m$与重复次数$N$关系")
plt.text(150, 6, '$m = clog_{2}N$')
plt.axis([0, 200, 0, 8])
plt.legend()
plt.show()
