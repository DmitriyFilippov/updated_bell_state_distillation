from definitions import *

def iceberg_simulation(input_error, n):

    stabilisers = ["X"*n, "Z"*n]

    out = kron(*[id for _ in range(n-2)], zerovector.T, plusvector.T)
    for i in range(2):
        out = np.matmul(mqpm(n, stabilisers[i], 0), out) * np.sqrt(2)
    
    out = np.matmul(kron(*[pauli_error_dict[x] for x in input_error]), out)

    for i in range(2):
        out = np.matmul(mqpm(n, stabilisers[i], 0), out)
    
    pout = np.trace(np.matmul(out, out.conjugate().T)) / 2**(n-2)
    out = np.matmul(kron(*[id for _ in range(n-2)], zerovector, plusvector), out) * (np.sqrt(2) ** 2)


    pid = np.abs(np.trace(conjugate(kron(*[id]*(n-2), out), duplicate_bell_state(bell, n-2)))) ** 2

    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0
    if pout > 1e-5 and pid < pout * 0.99999:
        for output_error in itertools.product(["I", "X", "Y", "Z"], repeat=n-2):
            if output_error.count("I") == n-2:
                continue
            perror = np.abs(np.trace(conjugate(kron(*[pauli_error_dict[x] for x in output_error], out), duplicate_bell_state(bell, n-2)))) ** 2
            if perror > 1e-5:
                totalpxerror += perror * output_error.count("X")
                totalpyerror += perror * output_error.count("Y")
                totalpzerror += perror * output_error.count("Z")
                totalperror += perror

    return pout, totalperror, totalpxerror, totalpyerror, totalpzerror, pid

print(iceberg_simulation(("I", "I", "I", "I", "I", "X"), 6))