"""
Author: Justin
Date: 2024-08-14 08:58:03
fielname: 
version: 
Description: 
LastEditTime: 2024-08-14 08:58:04
"""

from field import FieldElement
from polynomial import interpolate_poly, X, Polynomial
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


# 获取下一层的定义域
# 就是获取这一层的前一半的元素，然后再平方
def next_fri_domain(fri_domain):
    tmp = []
    # 获取fri_domain的前一半的元素，再平方
    for i in range(len(fri_domain) // 2):
        tmp.append(fri_domain[i] ** 2)
    return tmp


# 获取下一层的多项式
# 就是把多项式的奇数系数和偶数系数分别取出来，然后奇数系数乘以beta加上偶数系数
# 实现了多项式的降阶
def next_fri_polynomial(poly, beta):
    # 这里取出来了两个系数列表，一个是奇数项，一个是偶数项
    odd_coefficients = poly.poly[1::2]  # No need to fix this line.
    even_coefficients = poly.poly[::2]  # No need to fix this line either.
    # 系数列表减半了，所以多项式的阶也就减半了，实现了多项式的降阶
    # 然后在这里构造了新的多项式
    odd = beta * Polynomial(odd_coefficients)
    even = Polynomial(even_coefficients)
    return odd + even


# 获取FRI的下一层
def next_fri_layer(poly, domain, beta):
    next_poly = next_fri_polynomial(poly, beta)
    next_domain = next_fri_domain(domain)
    next_layer = []
    for x in next_domain:
        next_layer.append(next_poly(x))
    return next_poly, next_domain, next_layer


# 构建FRI Commitment
def FriCommit(cp, domain, cp_eval, cp_merkle, channel):
    fri_polys = [cp]
    fri_domains = [domain]
    fri_layers = [cp_eval]
    fri_merkles = [cp_merkle]
    while (
        fri_polys[-1].degree() > 0
    ):  # Replace this with the correct halting condition.
        beta = (
            channel.receive_random_field_element()
        )  # Change to obtain a random element from the channel.
        next_poly, next_domain, next_layer = next_fri_layer(
            fri_polys[-1], fri_domains[-1], beta
        )  # Fix to obtain the next FRI polynomial, domain, and layer.

        fri_polys.append(next_poly)
        fri_domains.append(next_domain)
        fri_layers.append(next_layer)

        fri_merkle = MerkleTree(next_layer)
        channel.send(fri_merkle.root)  # Fix to send the correct commitment.
        fri_merkles.append(fri_merkle)  # Fix to construct the correct Merkle tree.
    channel.send(str(fri_polys[-1].poly[0]))  # Fix to send the last layer's free term.
    return fri_polys, fri_domains, fri_layers, fri_merkles


"""
==================================================================
==============================Part.4==============================
==================================================================
"""
