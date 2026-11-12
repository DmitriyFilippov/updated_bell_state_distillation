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
