from definitions import *

def apply_storage_errors(state, n, p, xerrorvolume, yerrorvolume, zerrorvolume, count, target):
    if count < target or count >= target + (n * 3):
        return (state, p)
    i = count - target
    qubit = i // 3
    if i % 3 == 0:
        errorstring = ["I"] * n
        errorstring[qubit] = "X"
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, xerrorvolume)
    elif i % 3 == 1:
        errorstring = ["I"] * n
        errorstring[qubit] = "Y"
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, yerrorvolume)
    elif i % 3 == 2:
        errorstring = ["I"] * n
        errorstring[qubit] = "Z"
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, zerrorvolume)



def apply_mqpm_error(state, n, p, errorstring, errorvolume, count, target):
    if count == target:
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, errorvolume)
    else:
        return (state, p)


def faulty_mqpm(n, pauli_array, result, p, errorvolume, count, target):
    if count == target:
        state = mqpm(n, pauli_array, 1-result)
        return (state, errorvolume)
    else:
        state = mqpm(n, pauli_array, result)
        return (state, p)

def apply_x_error(state, n, p, qubit, errorvolume, count, target):
    if count == target:
        errorstring = ["I"] * n
        errorstring[qubit] = "X"
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, errorvolume)
    else:
        return (state, p)

def apply_z_error(state, n, p, qubit, errorvolume, count, target):
    if count == target:
        errorstring = ["I"] * n
        errorstring[qubit] = "Z"
        state = np.matmul(kron(*[pauli_error_dict[x] for x in errorstring]), state)
        return (state, errorvolume)
    else:
        return (state, p)