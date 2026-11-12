import itertools
import numpy as np
from typing import List 
import mpmath
from string import ascii_letters
from mpmath import mp

rcp_sqrt2 = 1 / np.sqrt(2)

def digits_to_int(digits: List[int]) -> int:
    """
    Converts a bitstring (list of 0s and 1s) to an integer
    """
    return int("".join([str(i) for i in digits]),2)

def conjugate(m1, m2):
    """
    Returns m2 * m1 * m2.H
    """
    return np.matmul(m2, np.matmul(m1, m2.conjugate().T))

def cnot(control, target, dim):
    """
    Produce matrix for a CNOT acting in a space with dimention dim with given control and target qubits
    """
    mat = np.zeros((2**dim, 2**dim), dtype=np.cdouble)
    for digits in itertools.product([0,1], repeat=dim):
        digits = list(digits)
        output_digits = digits.copy()
        if digits[control] == 1:
            output_digits[target] = 1-output_digits[target]
        input = digits_to_int(digits)
        output = digits_to_int(output_digits)
        mat[output, input] = 1
    return mat

def cz(control, target, dim):
    """
    Produce matrix for a CZ acting in a space with dimention dim with given control and target qubits
    """
    mat = np.zeros((2**dim, 2**dim), dtype=np.cdouble)
    for digits in itertools.product([0,1], repeat=dim):
        digits = list(digits)
        row = digits_to_int(digits)
        if digits[control] == 1 and digits[target] == 1:
            mat[row, row] = -1
        else:
            mat[row, row] = 1
    return mat

def swap(target, dim):
    """
    Produces a permutation matrix, where target indicated the position into which corresponding qubit should be swapped
    """
    mat = np.zeros((2**dim, 2**dim), dtype=np.cdouble)
    for digits in itertools.product([0,1], repeat=dim):
        digits = list(digits)
        output_digits = digits.copy()
        for i in range(dim):
            output_digits[target[i]] = digits[i]
        input = digits_to_int(digits)
        output = digits_to_int(output_digits)
        mat[output, input] = 1
    return mat

def kron(*args: np.matrix) -> np.matrix:
    """
    Produces Kroneker product of arguments
    """
    if (len(args) == 1):
        return args[0]
    out = np.kron(args[0], args[1])
    for j in range(2, len(args)):
        out = np.kron(out, args[j])
    return out

def duplicate_bell_state(bell, rep):
    """
    Duplicates a noisy bell state rep times
    """
    out = kron(*[bell]*rep)
    swaps = []
    for i in range(rep):
        swaps += [i, i+rep]
    out = conjugate(out, swap(swaps, rep*2))
    return out

# def mqpm(n, pauli_array, result):
#     paulis = []
#     paulis_t = []
#     identity_positions = []
#     for i in range(n):
#         pauli = pauli_array[i]
#         if pauli == "X":
#             paulis.append(had)
#             paulis_t.append(had)
#         elif pauli == "Y":
#             paulis.append(v)
#             paulis_t.append(v_t)
#         elif pauli == "Z":
#             paulis.append(id)
#             paulis_t.append(id)
#         elif pauli == "I":
#             paulis.append(id)
#             paulis_t.append(id)
#             identity_positions.append(i)

#     identities = len(identity_positions)

#     copy = np.zeros((2**(n), 2**(n)), dtype=np.cdouble)
#     copy_pi = np.zeros((2**(n), 2**(n)), dtype=np.cdouble)

#     for string in itertools.product([0, 1], repeat=identities):
        
#         zeroes = [0] * n
#         for j in range(identities):
#             zeroes[identity_positions[j]] = string[j]
#         i = digits_to_int(zeroes)
#         copy[i,i] = 1
#         copy_pi[i,i] = 1

#         ones = [1] * n
#         for j in range(identities):
#             ones[identity_positions[j]] = string[j]
#         j = digits_to_int(ones)

#         copy[j,j] = 1
#         copy_pi[j,j] = -1

#     state = kron(*paulis)
    
#     if result == 0:
#         state = np.matmul(copy, state)
#     else:
#         state = np.matmul(copy_pi, state)

#     state = np.matmul(kron(*paulis_t), state)
#     return state

def mqpm(n, pauli_array, result):
    paulis = []
    for pauli in pauli_array:
        if pauli == "X":
            paulis.append(had)
        elif pauli == "Y":
            paulis.append(v)
        elif pauli == "Z":
            paulis.append(id)
        elif pauli == "I":
            paulis.append(id)
    paulis.append(id)
    paulis_t = []
    for pauli in pauli_array:
        if pauli == "X":
            paulis_t.append(had)
        elif pauli == "Y":
            paulis_t.append(v_t)
        elif pauli == "Z":
            paulis_t.append(id)
        elif pauli == "I":
            paulis_t.append(id)
    paulis_t.append(id)
    

    state = kron(*[id for _ in range(n)], zerovector.T)
    state = np.matmul(kron(*paulis), state)
    for i in range(n):
        if pauli_array[i] != "I":
            state = np.matmul(cnot(i, n, n + 1), state) 
    state = np.matmul(kron(*paulis_t), state)
    if result == 0:
        state = np.matmul(kron(*[id for _ in range(n)], zerovector), state)
    else:
        state = np.matmul(kron(*[id for _ in range(n)], onevector), state)
    return state

def partial_trace(matrix, dim, indices):
    matrix = np.array(matrix, dtype=np.cdouble)
    matrix = np.reshape(matrix, [1] + [2] * 2 * dim)
    indices = np.sort(indices)
    rho_dim = 2 * dim
    for i, target_index in enumerate(indices):
        target_index = target_index - i
        state_indices = ascii_letters[1 : rho_dim - 2 * i + 1]
        state_indices = list(state_indices)

        target_letter = state_indices[target_index]
        state_indices[target_index + dim - i] = target_letter
        state_indices = "".join(state_indices)

        einsum_indices = f"a{state_indices}"
        matrix = np.einsum(einsum_indices, matrix)

    number_wires_sub = dim - len(indices)
    reduced_density_matrix = np.reshape(
        matrix, (1, 2**number_wires_sub, 2**number_wires_sub)
    )
    return reduced_density_matrix[0]

def factorize(matrix, dim):
    matrices = []
    for i in range(dim):
        indices = [x for x in range(dim)]
        indices.remove(i)
        matrices += [partial_trace(matrix, dim, indices)]
    return matrices


# Basic 1-qubit unitaries
x = np.matrix([[0,1],[1,0]], dtype=np.cdouble)
v = np.matrix([[0.5+0.5j,0.5-0.5j],[0.5-0.5j,0.5+0.5j]], dtype=np.cdouble)
v_t = np.matrix([[0.5-0.5j,0.5+0.5j],[0.5+0.5j,0.5-0.5j]], dtype=np.cdouble)

z = np.matrix([[1,0],[0,-1]], dtype=np.cdouble)
s = np.matrix([[1,0],[0,0+1j]], dtype=np.cdouble)
s_t = np.matrix([[1,0],[0,-1j]], dtype=np.cdouble)

y = np.matrix([[0,-1j],[1j,0]], dtype=np.cdouble)

id = np.matrix([[1,0],[0,1]], dtype=np.cdouble)

had = rcp_sqrt2 * np.matrix([[1,1],[1,-1]], dtype=np.cdouble)


pauli_error_dict = { 
    "I": id,
    "X": x,
    "Y": y,
    "Z": z
}

# Basic states as vectors


zerovector = np.matrix([1,0], dtype=np.cdouble)
onevector = np.matrix([0,1], dtype=np.cdouble)
plusvector = np.matrix([rcp_sqrt2, rcp_sqrt2], dtype=np.cdouble)
minusvector = np.matrix([rcp_sqrt2, -rcp_sqrt2], dtype=np.cdouble)
rightvector = np.matrix([rcp_sqrt2, rcp_sqrt2 * 1j], dtype=np.cdouble)
leftvector = np.matrix([rcp_sqrt2, -rcp_sqrt2 * 1j], dtype=np.cdouble)

comp_basis_vector_dict = { 
    0: zerovector,
    1: onevector,
}


had_basis_vector_dict = {
    0: plusvector,
    1: minusvector,
}

# Basic states as density matrices

zero = np.matmul(zerovector.conjugate().T, zerovector)
one = np.matmul(onevector.conjugate().T, onevector)
plus = np.matmul(plusvector.conjugate().T, plusvector)
minus = np.matmul(minusvector.conjugate().T, minusvector)
bell = cnot(1,0,2)*kron(zero, plus)*cnot(1,0,2).T


comp_basis_dict = {
    0: zero,
    1: one,
}

had_basis_dict = {
    0: plus,
    1: minus,
}