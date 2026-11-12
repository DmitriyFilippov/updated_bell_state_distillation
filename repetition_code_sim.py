from definitions import *

def z_repetition_simulation(input_error):


    out = kron(id, plusvector.T)

    out = np.matmul(mqpm(2, "ZZ", 0), out) * np.sqrt(2)
    
    out = np.matmul(kron(*[pauli_error_dict[x] for x in input_error]), out)

    out = np.matmul(mqpm(2, "ZZ", 0), out)

    pout = np.trace(np.matmul(out, out.conjugate().T)) / 2
    out = np.matmul(kron(id, plusvector), out) * (np.sqrt(2))


    pid = np.abs(np.trace(conjugate(kron(id, out), bell))) ** 2
    perror = pout - pid
    pxerror = np.abs(np.trace(conjugate(kron(x, out), bell))) ** 2
    pyerror = np.abs(np.trace(conjugate(kron(y, out), bell))) ** 2
    pzerror = np.abs(np.trace(conjugate(kron(z, out), bell))) ** 2


    return pout, perror, pxerror, pyerror, pzerror, pid

def x_repetition_simulation(input_error):


    out = kron(id, zerovector.T)

    out = np.matmul(mqpm(2, "XX", 0), out) * np.sqrt(2)
    
    out = np.matmul(kron(*[pauli_error_dict[x] for x in input_error]), out)

    out = np.matmul(mqpm(2, "XX", 0), out)

    pout = np.trace(np.matmul(out, out.conjugate().T)) / 2
    out = np.matmul(kron(id, zerovector), out) * (np.sqrt(2))


    pid = np.abs(np.trace(conjugate(kron(id, out), bell))) ** 2
    perror = pout - pid
    pxerror = np.abs(np.trace(conjugate(kron(x, out), bell))) ** 2
    pyerror = np.abs(np.trace(conjugate(kron(y, out), bell))) ** 2
    pzerror = np.abs(np.trace(conjugate(kron(z, out), bell))) ** 2


    return pout, perror, pxerror, pyerror, pzerror, pid


simulator = x_repetition_simulation

preject = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
prejectx = 0
prejecty = 0
prejectz = 0
#probability of different types of errors given in terms of polynomial in P_X, P_Y, P_Z
xs = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
ys = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
zs = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
total_count = [0, 0, 0, 0, 0, 0]

# Modeling of input errors. 5 qubits, 3 types of errors (+ no error) = 4^5 = 1024 input errors
for error in itertools.product(["I", "X", "Y", "Z"], repeat=2):
    input_error_name = "".join(error)
    input_error_size = 2 - input_error_name.count("I")
    pout, perror, pxerror, pyerror, pzerror, pid  = simulator(error)
    # Number of errors that happen on the input
    preject[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += 1-pout
    total_count[input_error_size - 1] +=  perror
    xs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pxerror
    ys[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pyerror
    zs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pzerror

print(preject)
print(total_count)

for xerr in range(3):
    for yerr in range(3):
        for zerr in range(3):
            if preject[xerr][yerr][zerr] > 1e-8:
                print(f"Rejections given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {preject[xerr][yerr][zerr]}")

for xerr in range(3):
    for yerr in range(3):
        for zerr in range(3):
            if xs[xerr][yerr][zerr] > 1e-8:
                print(f"X error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {xs[xerr][yerr][zerr]}")

for xerr in range(3):
    for yerr in range(3):
        for zerr in range(3):
            if ys[xerr][yerr][zerr] > 1e-8:
                print(f"Y error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {ys[xerr][yerr][zerr]}")

for xerr in range(3):
    for yerr in range(3):
        for zerr in range(3):
            if zs[xerr][yerr][zerr] > 1e-8:
                print(f"Z error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {zs[xerr][yerr][zerr]}")
