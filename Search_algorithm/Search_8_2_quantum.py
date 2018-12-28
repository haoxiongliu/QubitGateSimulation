import matplotlib.pyplot as plt
from pylab import mpl

from clsQubit import *
from gatelib import *

mpl.rcParams['font.sans-serif'] = ['FangSong']  # 指定默认字体
mpl.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题

n = 4
M, N = 2, 8
R = int(math.floor(np.pi * np.sqrt(N / M) / 4))

q_vec = QubitVec(n)
init = np.zeros(16)
init[1] = 1
q_vec.init_state('typein', tuple(init))

init_super = product_all((('H', n, 0), ('H', n, 1), ('H', n, 2), ('H', n, 3)))
grover_unit = product_all((('X', n, 3, (1,), (0,)), ('H', n, 0), ('H', n, 1),
                           ('H', n, 2), ('X', n, 0), ('X', n, 1), ('X', n, 2),
                           ('H', n, 2), ('X', n, 2, (0, 1)), ('H', n, 2),
                           ('X', n, 0), ('X', n, 1), ('X', n, 2),
                           ('H', n, 0), ('H', n, 1), ('H', n, 2)))
out_gate = product_all((('H', n, 3),))

iter_max = 10
pos_applytimes = list(np.zeros(iter_max))
for iter_times in range(0, iter_max):
    cir_U = init_super
    for k in range(0, iter_times):
        cir_U = grover_unit.dot(cir_U)
    cir_U = out_gate.dot(cir_U)

    q_vec.init_state('typein', tuple(init))
    q_vec.in_gate(cir_U)
    q_vec.calc_pos()
    pos_applytimes[iter_times] = \
        round(q_vec.calc_pos_state((0, 1, 1, 1)) + q_vec.calc_pos_state((0, 1, 0, 1)), 4)

print(pos_applytimes)
theta = 2 * np.arccos(np.sqrt((N - M) / N))
x = np.linspace(0, iter_max - 0.5, 100)
y = (np.sin((x + 1 / 2) * theta)) ** 2
plt.plot(x, y, linestyle='--', label=u'理论曲线')
plt.plot([i for i in range(0, iter_max)], pos_applytimes, 'ro', label=u'模拟数据')
plt.annotate('R = {}, p = {}'.format(R, pos_applytimes[R]),
             (R, pos_applytimes[R]), arrowprops=dict(arrowstyle='->'))
plt.ylabel(u"测得对应状态的概率")
plt.xlabel(u"grover迭代次数")
plt.title(u"M = {}, N = {}， R = {}".format(M, N, R))
plt.legend()
plt.show()
