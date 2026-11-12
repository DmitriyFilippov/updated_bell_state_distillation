from definitions import *

def perfect_simulation(input_error):

    stabilisers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]

    out = kron(id, zerovector.T, zerovector.T, zerovector.T, zerovector.T)
    for i in range(4):
        out = np.matmul(mqpm(5, stabilisers[i], 0), out) * np.sqrt(2)
    
    out = np.matmul(kron(*[pauli_error_dict[x] for x in input_error]), out)

    for i in range(4):
        out = np.matmul(mqpm(5, stabilisers[i], 0), out)

    pout = np.trace(np.matmul(out, out.conjugate().T)) / 2
    out = np.matmul(kron(id, zerovector, zerovector, zerovector, zerovector), out) * (np.sqrt(2) ** 4)


    pid = np.abs(np.trace(conjugate(kron(id, out), bell))) ** 2
    perror = pout - pid
    pxerror = np.abs(np.trace(conjugate(kron(x, out), bell))) ** 2
    pyerror = np.abs(np.trace(conjugate(kron(y, out), bell))) ** 2
    pzerror = np.abs(np.trace(conjugate(kron(z, out), bell))) ** 2


    return pout, perror, pxerror, pyerror, pzerror, pid

print(perfect_simulation(("Z", "I", "X", "X", "I")))