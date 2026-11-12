from definitions import *
from clifford_error_enumeration_definitions import (
    apply_mqpm_error,
    apply_storage_errors,
    faulty_mqpm,
)

def apply_iceberg_x_mqpm_errors(state, n, p, count, target):
    # totalerrors = n + (2/n) - 1
    (state, p) = apply_mqpm_error(state, n, p, "XX" + "I" * (n - 2), 4 + 2/3, count, target)
    (state, p) = apply_mqpm_error(state, n, p, "I" * (n - 2) + "XX", 4 + 2/3, count, target + 1)
    for i in range(n):
            errorstring = ["I"] * n
            errorstring[i] = "X"
            (state, p) = apply_mqpm_error(state, n, p, "".join(errorstring), 1/3, count, target + i)
    for i in range(0, n//2 - 3):
        break_after_qubits = (i + 1) * 2
        errorstring = ["X"] * break_after_qubits + ["I"] * (n-break_after_qubits)
        (state, p) = apply_mqpm_error(state, n, p, "".join(errorstring), 3 + 2/3, count, target + n + i)
    for i in range(0, n//2 - 2):
        errorstring = ["I"] * n
        errorstring[i * 2] = "X"
        errorstring[i * 2 + 1] = "X"
        (state, p) = apply_mqpm_error(state, n, p, "".join(errorstring), 2/3, count, target + n + i)

def apply_iceberg_z_mqpm_errors(state, n, p, count, target):
    # totalerrors = 2n - 3
    (state, p) = apply_mqpm_error(state, n, p, "Z" + "I" * (n - 1), 2 + 1/3, count, target)
    (state, p) = apply_mqpm_error(state, n, p, "I" * (n - 1) + "Z", 2 + 1/3, count, target + 1)
    for i in range(0, n-2):
            qubit = i + 1
            errorstring = ["I"] * n
            errorstring[qubit] = "Z"
            (state, p) = apply_mqpm_error(state, n, p, "".join(errorstring), 1/3, count, target + 2 + i)
    for i in range(0, n - 3):
        break_after_qubits = i + 2
        errorstring = ["Z"] * break_after_qubits + ["I"] * (n-break_after_qubits)
        (state, p) = apply_mqpm_error(state, n, p, "".join(errorstring), 1 + 2/3, count, target + n + i)

    return state, p

def iceberg_simulation(count, n):
    p=0

    state = kron(*[id for _ in range(n-2)], zerovector.T, plusvector.T)
    (state, p) = apply_storage_errors(state, n, p, 1/3, 1/3, 1/3, count, 0) # pyright: ignore[reportGeneralTypeIssues]
    


    (op,p) = faulty_mqpm(n, "X"*n, 0, p, 3*n-4, count, 15)
    state = np.matmul(op, state) * np.sqrt(2)

    
    

    (op,p) = faulty_mqpm(n, "Z"*n, 0, p, 2*n-1, count, 15)
    state = np.matmul(op, state) * np.sqrt(2)

    (state, p) = apply_storage_errors(state, n, p, 1/3, 1/3, 1/3, count, 27)  # pyright: ignore[reportGeneralTypeIssues]


    (op,p) = faulty_mqpm(n, "X"*n, 0, p, 3*n-4, count, 15)
    state = np.matmul(op, state) * np.sqrt(2)
    

    (op,p) = faulty_mqpm(n, "Z"*n, 0, p, 2*n-1, count, 15)
    state = np.matmul(op, state) * np.sqrt(2)


    (state, p) = apply_storage_errors(state, n, p, 1/3, 1/3, 1/3, count, 27)  # pyright: ignore[reportGeneralTypeIssues]


    pout = np.trace(np.matmul(state, state.conjugate().T)) / 2**(n-2)
    state = np.matmul(kron(id, zerovector, zerovector, zerovector, zerovector), state) * (np.sqrt(2) ** 4)

    pid = np.abs(np.trace(conjugate(kron(*[id]*(n-2), state), duplicate_bell_state(bell, n-2)))) ** 2


    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0
    if pout > 1e-5 and pid < pout * 0.99999:
        for stateput_error in itertools.product(["I", "X", "Y", "Z"], repeat=n-2):
            if stateput_error.count("I") == n-2:
                continue
            perror = np.abs(np.trace(conjugate(kron(*[pauli_error_dict[x] for x in stateput_error], state), duplicate_bell_state(bell, n-2)))) ** 2
            if perror > 1e-5:
                totalpxerror += perror * stateput_error.count("X")
                totalpyerror += perror * stateput_error.count("Y")
                totalpzerror += perror * stateput_error.count("Z")
                totalperror += perror

    return pout, totalperror, totalpxerror, totalpyerror, totalpzerror, pid, p
