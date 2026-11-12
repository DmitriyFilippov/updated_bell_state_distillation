from definitions import *
from clifford_error_enumeration_definitions import (
    apply_mqpm_error,
    apply_storage_errors,
    faulty_mqpm,
)

def z_repetition_simulation(count):
    p = 0

    state = kron(id, plusvector.T)
    (state, p) = apply_storage_errors(state, 2, p, 1/3, 1/3, 1/3, count, 0) # pyright: ignore[reportGeneralTypeIssues]


    (op,p) = faulty_mqpm(2, "ZZ", 0, p, 1, count, 6)
    state = np.matmul(op, state) * np.sqrt(2)
    (state, p) = apply_mqpm_error(state, 2, p, "ZI", 1, count, 7)
    
    (state, p) = apply_storage_errors(state, 2, p, 2/3, 2/3, 2/3, count, 8) # pyright: ignore[reportGeneralTypeIssues]

    (op,p) = faulty_mqpm(2, "ZZ", 0, p, 1, count, 14)
    state = np.matmul(op, state)
    (state, p) = apply_mqpm_error(state, 2, p, "ZI", 1, count, 15)

    (state, p) = apply_storage_errors(state, 2, p, 1/3, 1/3, 1/3, count, 16) # pyright: ignore[reportGeneralTypeIssues]

    pout = np.trace(np.matmul(state, state.conjugate().T)) / 2
    state = np.matmul(kron(id, plusvector), state) * (np.sqrt(2))


    pid = np.abs(np.trace(conjugate(kron(id, state), bell))) ** 2
    perror = pout - pid
    pxerror = np.abs(np.trace(conjugate(kron(x, state), bell))) ** 2
    pyerror = np.abs(np.trace(conjugate(kron(y, state), bell))) ** 2
    pzerror = np.abs(np.trace(conjugate(kron(z, state), bell))) ** 2


    return pout, perror, pxerror, pyerror, pzerror, pid, p

def x_repetition_simulation(count):
    p = 0

    state = kron(id, zerovector.T)

    (state, p) = apply_storage_errors(state, 2, p, 1/3, 1/3, 1/3, count, 0) # pyright: ignore[reportGeneralTypeIssues]
    
    
    (op,p) = faulty_mqpm(2, "XX", 0, p, 1, count, 6)
    state = np.matmul(op, state) * np.sqrt(2)
    (state, p) = apply_mqpm_error(state, 2, p, "XI", 1, count, 7)
    
    (state, p) = apply_storage_errors(state, 2, p, 2/3, 2/3, 2/3, count, 8) # pyright: ignore[reportGeneralTypeIssues]

    (op,p) = faulty_mqpm(2, "XX", 0, p, 1, count, 14)
    state = np.matmul(op, state)
    (state, p) = apply_mqpm_error(state, 2, p, "XI", 1, count, 15)

    (state, p) = apply_storage_errors(state, 2, p, 1/3, 1/3, 1/3, count, 16) # pyright: ignore[reportGeneralTypeIssues]

    pout = np.trace(np.matmul(state, state.conjugate().T)) / 2
    state = np.matmul(kron(id, zerovector), state) * (np.sqrt(2))


    pid = np.abs(np.trace(conjugate(kron(id, state), bell))) ** 2
    perror = pout - pid
    pxerror = np.abs(np.trace(conjugate(kron(x, state), bell))) ** 2
    pyerror = np.abs(np.trace(conjugate(kron(y, state), bell))) ** 2
    pzerror = np.abs(np.trace(conjugate(kron(z, state), bell))) ** 2


    return pout, perror, pxerror, pyerror, pzerror, pid, p


totalpreject = 0
totalperror = 0
totalpxerror = 0
totalpyerror = 0
totalpzerror = 0
totalpid = 0

for counter in range(22):
    pout, perror, pxerror, pyerror, pzerror, pid, p = z_repetition_simulation(counter)

    totalpreject += (1-pout) * p
    totalperror += perror * p
    totalpxerror += pxerror * p
    totalpyerror += pyerror * p
    totalpzerror += pzerror * p
    totalpid += pid * p

print(totalpreject, totalperror, totalpxerror, totalpyerror, totalpzerror, totalpid)