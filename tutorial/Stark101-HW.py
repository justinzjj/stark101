"""
Author: Justin
Date: 2024-08-14 08:58:03
fielname: 
version: 
Description: 
LastEditTime: 2024-08-14 08:58:04
"""

from field import FieldElement
from polynomial import interpolate_poly, X
from merkle import MerkleTree
from channel import Channel

"""
==================================================================
==============================Part.1==============================
==================================================================
"""

# 1. 使用群元素构建数列
a = [FieldElement(1), FieldElement(3141592)]
while len(a) < 1023:
    tmp = a[-1] * a[-1] + a[-2] * a[-2]
    a.append(tmp)

# 2. 使用生成元素构建数列
g = FieldElement.generator() ** (3 * 2**30 / 1024)
G = []
for i in range(1024):
    G.append(g**i)

# 3. 构建多项式
f = interpolate_poly(G[:-1], a)

# 4. 在Coset上进行评估
w = FieldElement.generator()
h = w ** ((3 * (2**30)) // 8192)
H = [h**i for i in range(8192)]
eval_domain = [w * x for x in H]
f_eval = [f(d) for d in eval_domain]

# 5. 构建Merkle树，构建commitment
f_merkle = MerkleTree(f_eval)

# 6. 使用Channel发送根节点
channel = Channel()
channel.send(f_merkle.root)


"""
==================================================================
==============================Part.2==============================
==================================================================
"""

# 1. 第一个约束
numer0 = f - 1
denom0 = X - 1
p0 = numer0 / denom0

# 2. 第二个约束
numer1 = f - 2338775057
denom1 = X - g**1022
p1 = numer1 / denom1

# 3. 第三个约束
# Construct a list `lst` of the linear terms (x-g**i):


numer2 = f(g**2 * X) - f(g * X) ** 2 - f**2
denom2 = (X**1024 - 1) / ((X - g**1021) * (X - g**1022) * (X - g**1023))
p2 = numer2 / denom2

print("deg p0 =", p0.degree())
print("deg p1 =", p1.degree())
print("deg p2 =", p2.degree())


# 4. 组合多项式
def get_CP(channel):
    alpha0 = channel.receive_random_field_element()
    alpha1 = channel.receive_random_field_element()
    alpha2 = channel.receive_random_field_element()
    return alpha0 * p0 + alpha1 * p1 + alpha2 * p2


# 5. 对组合多项式进行承诺
def CP_eval(channel):
    CP = get_CP(channel)
    result = []
    for d in eval_domain:
        result.append(CP(d))
    return result


# 实际上就是用Merkle树对CP进行承诺
# CP又是p0, p1, p2的组合多项式，把CP在陪集的定义域上进行赋值
channel = Channel()
CP_merkle = MerkleTree(CP_eval(channel))
channel.send(CP_merkle.root)

"""
==================================================================
==============================Part.3==============================
==================================================================
"""
