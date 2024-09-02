"""
Author: Justin
Date: 2024-09-02 09:36:36
fielname: 
version: 
Description: 
LastEditTime: 2024-09-02 09:37:38
"""

from merkle import verify_decommitment
from field import FieldElement
from channel import Channel


# 返回下一层的域，也就是长度折半，但是每个元素都是原来的平方
def next_fri_domain(domain):
    return [x**2 for x in domain[: len(domain) // 2]]


class DecommitmentData:
    # 这里构造解密的数据结构，包括：
    # 索引
    # 值
    # Merkle树的路径
    # Merkle树的根
    def __init__(self, index, value, authentication_path, merkle_root):
        self.index = index
        self.value = value
        self.authentication_path = authentication_path
        self.merkle_root = merkle_root


class Proof:
    # 这里构造证明的数据结构，包括：
    # 前三项式Trace的证明
    x_proof: DecommitmentData
    gx_proof: DecommitmentData
    g2x_proof: DecommitmentData
    cp_proof: list[DecommitmentData]  # CP解密的证明
    final_values: list[FieldElement]  # 最终的那一个常数
    trace_domain: list[FieldElement]  # Trace的域
    lde_domain: list[FieldElement]  # LDE的域

    def __init__(
        self,
        x_proof: DecommitmentData,
        gx_proof: DecommitmentData,
        g2x_proof: DecommitmentData,
        cp_proof: list[DecommitmentData],
        final_values: list[FieldElement],
        trace_domain: list[FieldElement],
        lde_domain: list[FieldElement],
    ):
        self.x_proof = x_proof
        self.gx_proof = gx_proof
        self.g2x_proof = g2x_proof
        self.cp_proof = cp_proof if cp_proof is not None else []
        self.final_values = final_values
        self.trace_domain = trace_domain
        self.lde_domain = lde_domain

    # 1. 校验merkle证明
    def check_merkle_proof(self, merkle_proof_valid, lde_tree, fri_merkles):
        self.check_merkle_proof_Trace(merkle_proof_valid, lde_tree)
        self.check_merkle_proof_CP(merkle_proof_valid, fri_merkles)

    # 验证前三项 Trace的证明
    def check_merkle_proof_Trace(self, merkle_proof_valid, lde_tree):
        # 校验值 并且校验路径
        merkle_proof_valid.append(
            verify_decommitment(
                self.x_proof.index,
                self.x_proof.value,
                self.x_proof.authentication_path,
                self.x_proof.merkle_root,
            )
            & (self.x_proof.merkle_root == lde_tree.root)
        )
        merkle_proof_valid.append(
            verify_decommitment(
                self.gx_proof.index,
                self.gx_proof.value,
                self.gx_proof.authentication_path,
                self.gx_proof.merkle_root,
            )
            & (self.gx_proof.merkle_root == lde_tree.root)
        )
        merkle_proof_valid.append(
            verify_decommitment(
                self.g2x_proof.index,
                self.g2x_proof.value,
                self.g2x_proof.authentication_path,
                self.g2x_proof.merkle_root,
            )
            & (self.g2x_proof.merkle_root == lde_tree.root)
        )

    # 校验各层的CP证明和路径
    def check_merkle_proof_CP(self, merkle_proof_valid, fri_merkles):
        for i in range(len(self.cp_proof) // 2):
            proof = self.cp_proof[i * 2]
            merkle_proof_valid.append(
                verify_decommitment(
                    proof.index,
                    proof.value,
                    proof.authentication_path,
                    proof.merkle_root,
                )
                & (proof.merkle_root == fri_merkles[i].root)
            )
            proof = self.cp_proof[i * 2 + 1]
            merkle_proof_valid.append(
                verify_decommitment(
                    proof.index,
                    proof.value,
                    proof.authentication_path,
                    proof.merkle_root,
                )
                & (proof.merkle_root == fri_merkles[i].root)
            )

    # 3. 检查值之间的关系
    def check_relationship(self, alphas, channel):
        # check cp0(x), f(x), f(g*x), f(g^2*x) relationship
        x_val = self.lde_domain[self.x_proof.index]  # x = lde_domain[index]
        p0_val = (self.x_proof.value - 1) / (x_val - 1)  # p0 = (f(x) - 1)/(x-g^0)
        p1_val = (self.x_proof.value - 2338775057) / (
            x_val - self.trace_domain[1022]
        )  # p1 = (f(x) - 2338775057 )/(x-g^1022)

        p2_denominator = (x_val**1024 - 1) / (
            (x_val - self.trace_domain[1021])
            * (x_val - self.trace_domain[1022])
            * (x_val - self.trace_domain[1023])
        )
        p2_val = (
            self.g2x_proof.value - self.gx_proof.value**2 - self.x_proof.value**2
        ) / p2_denominator  # p2 = (f(g^2*x) - f(gx)^2 - f(x)^2)/[(x^1024 - 1)/((x-g^1021)*(x-g^1022)*(x-g^1023))]

        calculated_cp_val = p0_val * alphas[0] + p1_val * alphas[1] + p2_val * alphas[2]
        trace_cp_valid = calculated_cp_val == self.cp_proof[0].value

        # check cpi(x), cpi(-x), cp{i+1}(x^2) relationship
        domain = self.lde_domain
        cp_layers_valid = []
        for i in range(len(self.cp_proof) // 2):
            if i < len(self.cp_proof) // 2 - 1:
                next_cp_value = self.cp_proof[(i + 1) * 2].value  #
            else:
                next_cp_value = self.final_values[0]

            channel.send(self.cp_proof[i * 2].merkle_root)
            beta = channel.receive_random_field_element()
            idx = self.cp_proof[i * 2].index
            x_val = domain[idx]
            g_x_square_val = (
                self.cp_proof[i * 2].value + self.cp_proof[i * 2 + 1].value
            ) / 2
            h_x_square_val = (
                self.cp_proof[i * 2].value - self.cp_proof[i * 2 + 1].value
            ) / (2 * x_val)
            calculated_next_cp_eval = g_x_square_val + beta * h_x_square_val
            cp_layers_valid.append(calculated_next_cp_eval == next_cp_value)

            # update domain and next_cp_value
            domain = next_fri_domain(domain)
        return trace_cp_valid, cp_layers_valid

    # 验证
    def verify(self, final_value: FieldElement, lde_tree, fri_merkles):
        # 在整个验证过程中，需要验证以下内容：
        # 1. 所有的Merkle证明都是有效的
        # 2. 所有各自的根是相同的
        # 3. 值之间的关系
        # 4. 兄弟节点之间的路径正确

        # 1. 校验merkle证明
        merkle_proof_valid = []
        self.check_merkle_proof(merkle_proof_valid, lde_tree, fri_merkles)

        # 2. 校验各自最终的根是相同的：LDE、CP、Final Value
        same_lde_root = (
            self.x_proof.merkle_root
            == self.gx_proof.merkle_root
            == self.g2x_proof.merkle_root
        )
        same_cp_root = []
        for i in range(len(self.cp_proof) // 2):
            same_cp_root.append(
                self.cp_proof[i * 2].merkle_root == self.cp_proof[i * 2 + 1].merkle_root
            )
        same_final_values = []
        for i in range(len(self.final_values)):
            same_final_values.append(final_value == self.final_values[i])

        # 3. 值之间的关系
        channel = Channel()
        channel.send(self.x_proof.merkle_root)
        alphas = channel.derive_alphas(3)
        trace_cp_valid, cp_layers_valid = self.check_relationship(alphas, channel)

        # 4. 兄弟节点之间的路径正确
        # check sibling cp has same merkle root and lde_domain[idx]^2 == lde_domain[sibling_idx]^2
        same_sibling_square = []
        domain = self.lde_domain
        for i in range(len(self.cp_proof) // 2):
            idx = self.cp_proof[i * 2].index
            sibling_idx = self.cp_proof[i * 2 + 1].index
            same_sibling_square.append(domain[idx] ** 2 == domain[sibling_idx] ** 2)

            domain = next_fri_domain(domain)

        # 最终校验一下
        valid = (
            all(merkle_proof_valid)
            & same_lde_root
            & all(same_cp_root)
            & all(same_final_values)
            & trace_cp_valid
            & all(cp_layers_valid)
            & all(same_sibling_square)
        )
        return valid
