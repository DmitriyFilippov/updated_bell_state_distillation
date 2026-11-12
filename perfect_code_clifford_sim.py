from definitions import *
from clifford_error_enumeration_definitions import (
    apply_mqpm_error,
    apply_storage_errors,
    faulty_mqpm,
)

def perfect_clifford_error_simulation(count):
    p = 0 # probability of error factor

    state = kron(id, zerovector.T, zerovector.T, zerovector.T, zerovector.T)
    (state, p) = apply_storage_errors(state, 5, p, 1/3, 1/3, 1/3, count, 0) # pyright: ignore[reportGeneralTypeIssues]
    

    (op,p) = faulty_mqpm(5, "XZZXI", 0, p, 6, count, 15)
    state = np.matmul(op, state) * np.sqrt(2)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "XIIII", 2 + 1/3, count, 16)
    (state, p) = apply_mqpm_error(state, 5, p, "IZIII", 1 + 1/3, count, 17)
    (state, p) = apply_mqpm_error(state, 5, p, "IIZII", 1 + 1/3, count, 18)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIXI", 2/3, count, 19)
    (state, p) = apply_mqpm_error(state, 5, p, "XIIXI", 1/3, count, 20)



    (op,p) = faulty_mqpm(5, "ZXIXZ", 0, p, 6, count, 21)
    state = np.matmul(op, state) * np.sqrt(2)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IXIII", 2 + 1/3, count, 22)
    (state, p) = apply_mqpm_error(state, 5, p, "ZIIII", 1 + 1/3, count, 23)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIZ", 1 + 1/3, count, 24)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIXI", 1/3, count, 25)
    (state, p) = apply_mqpm_error(state, 5, p, "IXIXI", 1/3, count, 26)



    #storage errors
    (state, p) = apply_storage_errors(state, 5, p, 1/3, 1/3, 1/3, count, 27) # pyright: ignore[reportGeneralTypeIssues]



    (op,p) = faulty_mqpm(5, "IXZZX", 0, p, 6, count, 27)
    state = np.matmul(op, state) * np.sqrt(2)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IIZII", 2 + 1/3, count, 43)
    (state, p) = apply_mqpm_error(state, 5, p, "IXIII", 1 + 1/3, count, 44)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIX", 1 + 1/3, count, 45)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZI", 1/3, count, 46)
    (state, p) = apply_mqpm_error(state, 5, p, "IIZZI", 1/3, count, 47)




    (op,p) = faulty_mqpm(5, "XIXZZ", 0, p, 6, count, 48)
    state = np.matmul(op, state) * np.sqrt(2)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIZ", 2 + 1/3, count, 49)
    (state, p) = apply_mqpm_error(state, 5, p, "XIIII", 1 + 1/3, count, 50)
    (state, p) = apply_mqpm_error(state, 5, p, "IIXII", 1 + 1/3, count, 51)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZI", 1/3, count, 52)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZZ", 1/3, count, 53)




    #storage errors
    (state, p) = apply_storage_errors(state, 5, p, 2/3, 2/3, 2/3, count, 54) # pyright: ignore[reportGeneralTypeIssues]



    (op,p) = faulty_mqpm(5, "XZZXI", 0, p, 6, count, 69)
    state = np.matmul(op, state)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "XIIII", 2 + 1/3, count, 70)
    (state, p) = apply_mqpm_error(state, 5, p, "IZIII", 1 + 1/3, count, 71)
    (state, p) = apply_mqpm_error(state, 5, p, "IIZII", 1 + 1/3, count, 72)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIXI", 1/3, count, 73)
    (state, p) = apply_mqpm_error(state, 5, p, "XIIXI", 1/3, count, 74)



    (op,p) = faulty_mqpm(5, "ZXIXZ", 0, p, 6, count, 75)
    state = np.matmul(op, state)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IXIII", 2 + 1/3, count, 76)
    (state, p) = apply_mqpm_error(state, 5, p, "ZIIII", 1 + 1/3, count, 77)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIZ", 1 + 1/3, count, 78)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIXI", 1/3, count, 79)
    (state, p) = apply_mqpm_error(state, 5, p, "IXIXI", 1/3, count, 80)


    #storage errors
    (state, p) = apply_storage_errors(state, 5, p, 1/3, 1/3, 1/3, count, 81) # pyright: ignore[reportGeneralTypeIssues]



    (op,p) = faulty_mqpm(5, "IXZZX", 0, p, 6, count, 96)
    state = np.matmul(op, state)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IIZII", 2 + 1/3, count, 97)
    (state, p) = apply_mqpm_error(state, 5, p, "IXIII", 1 + 1/3, count, 98)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIX", 1 + 1/3, count, 99)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZI", 1/3, count, 100)
    (state, p) = apply_mqpm_error(state, 5, p, "IIZZI", 1/3, count, 101)



    (op,p) = faulty_mqpm(5, "XIXZZ", 0, p, 6, count, 102)
    state = np.matmul(op, state)
    #lattice surgery decoding errors
    (state, p) = apply_mqpm_error(state, 5, p, "IIIIZ", 2 + 1/3, count, 103)
    (state, p) = apply_mqpm_error(state, 5, p, "XIIII", 1 + 1/3, count, 104)
    (state, p) = apply_mqpm_error(state, 5, p, "IIXII", 1 + 1/3, count, 105)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZI", 1/3, count, 106)
    (state, p) = apply_mqpm_error(state, 5, p, "IIIZZ", 1/3, count, 107)


    #storage errors
    (state, p) = apply_storage_errors(state, 5, p, 1/3, 1/3, 1/3, count, 108) # pyright: ignore[reportGeneralTypeIssues]


    pout = np.trace(np.matmul(state, state.conjugate().T)) / 2
    state = np.matmul(kron(id, zerovector, zerovector, zerovector, zerovector), state) * (np.sqrt(2) ** 4)

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

for counter in range(123):
    pout, perror, pxerror, pyerror, pzerror, pid, p = perfect_clifford_error_simulation(counter)

    totalpreject += (1-pout) * p
    totalperror += perror * p
    totalpxerror += pxerror * p
    totalpyerror += pyerror * p
    totalpzerror += pzerror * p
    totalpid += pid * p

print(totalpreject, totalperror, totalpxerror, totalpyerror, totalpzerror, totalpid)