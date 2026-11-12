from definitions import *

def iceberg_simulation(error, n):

    stabilisers = ["X"*n, "Z"*n]

    totalpout = 0
    totalpid = 0
    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0

    for syndrome in itertools.product([0, 1], repeat=2):
        state = kron(*[id for _ in range(n-2)], zerovector.T, plusvector.T)
        for i in range(2):
            state = np.matmul(mqpm(n, stabilisers[i], syndrome[i]), state)
        
        state = np.matmul(kron(*[pauli_error_dict[x] for x in error]), state)

        for i in range(2):
            state = np.matmul(mqpm(n, stabilisers[i], syndrome[i]), state)
        
        
        for measurements in itertools.product([0, 1], repeat=2):
            out = np.matmul(kron(*[id for _ in range(n-2)], comp_basis_vector_dict[measurements[0]], had_basis_vector_dict[measurements[1]]), state)
    
            if measurements[0]:
                out = np.matmul(kron(*[x for _ in range(n-2)]), out)
            if measurements[1]:
                out = np.matmul(kron(*[z for _ in range(n-2)]), out)

            pout = np.trace(np.matmul(out, out.conjugate().T)) / (2**(n-2))
            if pout > 1e-5:
                for error2 in itertools.product(["I", "X", "Y", "Z"], repeat=n-2):
                    perror = np.abs(np.trace(conjugate(kron(*[pauli_error_dict[x] for x in error2], out), duplicate_bell_state(bell, n-2)))) ** 2
                    if perror > 1e-5:
                        print(f"Error: {error2}, Probability: {perror}")

for error in itertools.product(["I", "X", "Y", "Z"], repeat=6):
    print(f"Inpit Error: {error}")
    iceberg_simulation(error, 6)