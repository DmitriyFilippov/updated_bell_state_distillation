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



n=8

preject = [[[0 for _ in range(n+1)] for _ in range(n+1)] for _ in range(n+1)]
prejectx = 0
prejecty = 0
prejectz = 0
#probability of different types of errors given in terms of polynomial in P_X, P_Y, P_Z
xs = [[[0 for _ in range(n+1)] for _ in range(n+1)] for _ in range(n+1)]
ys = [[[0 for _ in range(n+1)] for _ in range(n+1)] for _ in range(n+1)]
zs = [[[0 for _ in range(n+1)] for _ in range(n+1)] for _ in range(n+1)]
total_count = [0, 0, 0, 0, 0, 0]

# Modeling of input errors. 5 qubits, 3 types of errors (+ no error) = 4^5 = 1024 input errors
for error in itertools.product(["I", "X", "Y", "Z"], repeat=n):
    input_error_name = "".join(error)
    input_error_size = n - input_error_name.count("I")
    if input_error_size > 2:
        continue
    pout, perror, pxerror, pyerror, pzerror, pid  = iceberg_simulation(error, n)
    # Number of errors that happen on the input
    preject[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += 1-pout
    total_count[input_error_size - 1] +=  perror
    xs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pxerror
    ys[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pyerror
    zs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pzerror

print(preject)
print(total_count)

for xerr in range(n+1):
    for yerr in range(n+1):
        for zerr in range(n+1):
            if preject[xerr][yerr][zerr] > 1e-8:
                print(f"Rejections given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {preject[xerr][yerr][zerr]}")

for xerr in range(n+1):
    for yerr in range(n+1):
        for zerr in range(n+1):
            if xs[xerr][yerr][zerr] > 1e-8:
                print(f"X error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {xs[xerr][yerr][zerr]}")

for xerr in range(n+1):
    for yerr in range(n+1):
        for zerr in range(n+1):
            if ys[xerr][yerr][zerr] > 1e-8:
                print(f"Y error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {ys[xerr][yerr][zerr]}")

for xerr in range(n+1):
    for yerr in range(n+1):
        for zerr in range(n+1):
            if zs[xerr][yerr][zerr] > 1e-8:
                print(f"Z error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {zs[xerr][yerr][zerr]}")
