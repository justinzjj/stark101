"""
Author: Justin
Date: 2024-08-19 10:54:07
fielname: 
version: 
Description: 
LastEditTime: 2024-08-19 10:54:08
"""

import time
from field import FieldElement
from polynomial import interpolate_poly, X, Polynomial
from merkle import MerkleTree
from channel import Channel


# 辅助生产标题
def print_section_header(section_name, total_width=66):
    # Calculate the number of equals signs needed on each side
    name_length = len(section_name)
    side_width = (
        total_width - name_length - 2
    ) // 2  # Subtract 2 for the spaces around the section name

    # Create the header lines
    top_bottom_line = "=" * total_width
    middle_line = f"{'=' * side_width} {section_name} {'=' * side_width}"

    # Adjust middle line if the total width is odd
    if len(middle_line) < total_width:
        middle_line += "="

    # Print the header
    print(f"\n{top_bottom_line}\n{middle_line}\n{top_bottom_line}\n")


def gengerating_trace(channel):
    # 1. 使用群元素构建数列
    a = [FieldElement(1), FieldElement(3141592)]
    while len(a) < 1023:
        tmp = a[-1] * a[-1] + a[-2] * a[-2]
        a.append(tmp)

    # 2. 使用生成元素构建群元素
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
    channel.send(f_merkle.root)
    return (g, G, eval_domain, f, f_eval, f_merkle)


def generating_composition_polynomial(g, eval_domain, f, channel):

    # 1. 第一个约束
    numer0 = f - 1
    denom0 = X - 1
    p0 = numer0 / denom0
    # 2. 第二个约束
    numer1 = f - 2338775057
    denom1 = X - g**1022
    p1 = numer1 / denom1
    # 3. 第三个约束
    numer2 = f(g**2 * X) - f(g * X) ** 2 - f**2
    denom2 = (X**1024 - 1) / ((X - g**1021) * (X - g**1022) * (X - g**1023))
    p2 = numer2 / denom2

    cp = get_cp(channel, p0, p1, p2)

    cp_eval_result = cp_eval(channel, eval_domain, cp)

    cp_merkle = MerkleTree(cp_eval_result)
    channel.send(cp_merkle.root)

    return cp, cp_merkle, cp_eval_result, channel


def cp_eval(channel, eval_domain, CP):
    result = []
    for d in eval_domain:
        result.append(CP(d))
    return result


# 4. 组合多项式
def get_cp(channel, p0, p1, p2):
    alpha0 = channel.receive_random_field_element()
    alpha1 = channel.receive_random_field_element()
    alpha2 = channel.receive_random_field_element()
    return alpha0 * p0 + alpha1 * p1 + alpha2 * p2


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
    return fri_layers, fri_merkles


# 获取FRI的下一层
def next_fri_layer(poly, domain, beta):
    next_poly = next_fri_polynomial(poly, beta)
    next_domain = next_fri_domain(domain)
    next_layer = []
    for x in next_domain:
        next_layer.append(next_poly(x))
    return next_poly, next_domain, next_layer


def next_fri_polynomial(poly, beta):
    # 这里取出来了两个系数列表，一个是奇数项，一个是偶数项
    odd_coefficients = poly.poly[1::2]  # No need to fix this line.
    even_coefficients = poly.poly[::2]  # No need to fix this line either.
    # 系数列表减半了，所以多项式的阶也就减半了，实现了多项式的降阶
    # 然后在这里构造了新的多项式
    odd = beta * Polynomial(odd_coefficients)
    even = Polynomial(even_coefficients)
    return odd + even


def next_fri_domain(fri_domain):
    tmp = []
    # 获取fri_domain的前一半的元素，再平方
    for i in range(len(fri_domain) // 2):
        tmp.append(fri_domain[i] ** 2)
    return tmp


# 解码fri
def decommit_fri(channel, f_eval, f_merkle, fri_layers, fri_merkles):
    for query in range(3):
        decommit_on_query(
            channel.receive_random_int(0, 8191 - 16),
            channel,
            f_eval,
            f_merkle,
            fri_layers,
            fri_merkles,
        )

    return channel


# 解码查询
def decommit_on_query(idx, channel, f_eval, f_merkle, fri_layers, fri_merkles):
    # 输入索引和channel
    channel.send(str(f_eval[idx]))  # 发送f(x)的值
    channel.send(str(f_merkle.get_authentication_path(idx)))  # 发送f(x)的验证路径
    channel.send(str(f_eval[idx + 8]))  # 发送f(gx)的值
    channel.send(str(f_merkle.get_authentication_path(idx + 8)))  # 发送f(gx)的验证路径
    channel.send(str(f_eval[idx + 16]))  # 发送f(g^2x)的值
    channel.send(
        str(f_merkle.get_authentication_path(idx + 16))
    )  # 发送f(g^2x)的验证路径
    decommit_on_fri_layers(
        idx,
        channel,
        fri_layers,
        fri_merkles,
    )


# 解码frilayers
def decommit_on_fri_layers(idx, channel, fri_layers, fri_merkles):
    # 传入channel，和idx 索引
    for layer, merkle in zip(fri_layers[:-1], fri_merkles[:-1]):
        length = len(layer)
        idx = idx % length  # 获取了cp（x）的索引
        sib_idx = (idx + length // 2) % length  # 获取了cp（-x）的索引
        channel.send(str(layer[idx]))
        channel.send(str(merkle.get_authentication_path(idx)))
        channel.send(str(layer[sib_idx]))
        channel.send(str(merkle.get_authentication_path(sib_idx)))
    channel.send(
        str(fri_layers[-1][0])
    )  # 最后一层的元素，这里只发送一个 因为最后一层只有一个元素


def main():
    start = time.time()

    channel = Channel()
    # Generating the trace...
    start_all = start
    print_section_header("Generating the trace...")
    g, G, eval_domain, f, f_eval, f_merkle = gengerating_trace(channel)
    print(f"{time.time() - start}s")

    # Generating the composition polynomial
    start = time.time()
    print_section_header("Generating the composition polynomial and the FRI layers")

    cp, cp_merkle, cp_eval_result, channel = generating_composition_polynomial(
        g, eval_domain, f, channel
    )

    fri_layers, fri_merkles = FriCommit(
        cp, eval_domain, cp_eval_result, cp_merkle, channel
    )

    print(f"{time.time() - start}s")

    # Generating queries and decommitments
    start = time.time()
    print_section_header("Generating queries and decommitments...")
    channel = decommit_fri(channel, f_eval, f_merkle, fri_layers, fri_merkles)
    print(f"{time.time() - start}s")

    start = time.time()
    print_section_header(channel.proof)
    print_section_header(f"Overall time: {time.time() - start_all}s")
    print_section_header(
        f"Uncompressed proof length in characters: {len(str(channel.proof))}"
    )


if __name__ == "__main__":
    main()
