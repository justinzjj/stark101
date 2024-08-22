"""
Author: Justin
Date: 2024-08-22 09:03:14
fielname: 
version: 
Description: 
LastEditTime: 2024-08-22 09:04:26
"""

import time
from field import FieldElement
from polynomial import interpolate_poly, X, Polynomial, prod
from merkle import MerkleTree
from channel import Channel


def get_cp(channel, p0, p1, p2):
    alpha0 = channel.receive_random_field_element()
    alpha1 = channel.receive_random_field_element()
    alpha2 = channel.receive_random_field_element()
    return alpha0 * p0 + alpha1 * p1 + alpha2 * p2


# 1. get C
c = [FieldElement(1)]
for i in range(128):
    c.append(c[i] * 2 * (2 * i + 1) / (i + 2))

print("c[100] = " + str(c[100]))
# c[100]=1555040254

# 2. get G, G size is 128
g = FieldElement.generator() ** (3 * 2**30 / 128)
G = []
for i in range(128):
    G.append(g**i)

# 3. get f
f = interpolate_poly(G, c[:-1])


# 4. get H, H size is 1024
w = FieldElement.generator()
h = w ** ((3 * (2**30)) // 1024)
H = [h**i for i in range(1024)]
eval_domain = [w * x for x in H]

# 5.
T = interpolate_poly([g**i for i in range(128)], [FieldElement(p) for p in range(128)])

# 6. constrait p0,p1,p2
numer0 = f - 1
denom0 = X - 1
p0 = numer0 / denom0
numer1 = f - 1555040254
denom1 = X - g**100
p1 = numer1 / denom1
numer2 = f(g * X) * (T(X) + 2) - (2 * (2 * T(X) + 1) * f(X))
denom2 = (X**128 - 1) / (X - g**127)
p2 = numer2 / denom2

channel = Channel()
cp = get_cp(channel, p0, p1, p2)

print("cp(h*888) = " + str(cp(h**888)))
print("cp(g*111) = " + str(cp(g**111)))
