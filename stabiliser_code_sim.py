from definitions import *

def stabiliser_error_simulation(n, k, stabilisers, error, decoder):
    totalpout = 0
    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0
    for syndrome in itertools.product([0, 1], repeat=len(stabilisers)):
        state = kron(*[id for _ in range(k)], *[zerovector.T for _ in range(n-k)])
        for i in range(len(stabilisers)):
            state = np.matmul(mqpm(5, stabilisers[i], syndrome[i]), state)
        
        state = np.matmul(kron(*[pauli_error_dict[x] for x in error]), state)

        for i in range(len(stabilisers)):
            state = np.matmul(mqpm(5, stabilisers[i], syndrome[i]), state)

        for measurements in itertools.product([0, 1], repeat=n-k):
            out = np.matmul(kron(*[id for _ in range(k)], *[comp_basis_vector_dict[x] for x in measurements]), state)

            out = decoder(out, measurements)

            pout = np.trace(np.matmul(out, out.conjugate().T)) / (2**(n-2))
            pid = np.abs(np.trace(conjugate(kron(id, out), bell))) ** 2
            pxerror = np.abs(np.trace(conjugate(kron(x, out), bell))) ** 2
            pyerror = np.abs(np.trace(conjugate(kron(y, out), bell))) ** 2
            pzerror = np.abs(np.trace(conjugate(kron(z, out), bell))) ** 2

            totalpout += pout
            totalperror += (pxerror + pyerror + pzerror)
            totalpxerror += pxerror
            totalpyerror += pyerror
            totalpzerror += pzerror

    return totalpout, totalperror, totalpxerror, totalpyerror, totalpzerror

def fivequbit_decoder(state, syndrome):
    out = state
    if sum(syndrome) % 2:
        out = np.matmul(x, out)
    if (syndrome[0] + syndrome[3]) % 2:
        out = np.matmul(z, out)
    return out


def fivequbit_simulation(error):

    stabilisers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]

    totalpout = 0
    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0
    
    for syndrome in itertools.product([0, 1], repeat=4):
        state = kron(id, zerovector.T, zerovector.T, zerovector.T, zerovector.T)
        for i in range(4):
            state = np.matmul(mqpm(5, stabilisers[i], syndrome[i]), state)
        
        state = np.matmul(kron(*[pauli_error_dict[x] for x in error]), state)

        for i in range(4):
            state = np.matmul(mqpm(5, stabilisers[i], syndrome[i]), state)
        
        
        for measurements in itertools.product([0, 1], repeat=4):
            out = np.matmul(kron(id, *[comp_basis_vector_dict[x] for x in measurements]), state)
    
            if sum(measurements) % 2:
                out = np.matmul(x, out)
            if (measurements[0] + measurements[3]) % 2:
                out = np.matmul(z, out)


            pid = np.abs(np.trace(conjugate(kron(id, out), bell))) ** 2
            pxerror = np.abs(np.trace(conjugate(kron(x, out), bell))) ** 2
            pyerror = np.abs(np.trace(conjugate(kron(y, out), bell))) ** 2
            pzerror = np.abs(np.trace(conjugate(kron(z, out), bell))) ** 2

            totalpout += (pid + pxerror + pyerror + pzerror)
            totalperror += (pxerror + pyerror + pzerror)
            totalpxerror += pxerror
            totalpyerror += pyerror
            totalpzerror += pzerror
            
            
    if totalpout <= 1e-12:
        totalpout = 0
    if totalperror <= 1e-12:
        totalperror = 0
    if totalpxerror <= 1e-12:
        totalpxerror = 0
    if totalpyerror <= 1e-12:
        totalpyerror = 0
    if totalpzerror <= 1e-12:
        totalpzerror = 0
    return totalpout, totalperror, totalpxerror, totalpyerror, totalpzerror


def fivequbit_twoway_simulation(error):

    stabilisers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    
    totalpout = 0
    totalperror = 0
    totalpxerror = 0
    totalpyerror = 0
    totalpzerror = 0
    for syndrome in itertools.product([0, 1], repeat=4):
        state = duplicate_bell_state(bell, 5)
        state = conjugate(state, kron(*[pauli_error_dict[x] for x in error], id, id, id, id, id))
        for i in range(4):
            state = conjugate(state, kron(mqpm(5, stabilisers[i], syndrome[i]), mqpm(5, stabilisers[i], syndrome[i])))
        
        if np.trace(state).real < 1e-12:
            continue

        
        for measurements in itertools.product([0, 1], repeat=8):
            out = conjugate(state, kron(id, *[comp_basis_vector_dict[x] for x in measurements[0:4]], id, *[comp_basis_vector_dict[x] for x in measurements[4:8]]))
        
            if sum(measurements) % 2:
                out = conjugate(out, kron(x, id))
            if (measurements[0] + measurements[3] + measurements[4] + measurements[7]) % 2:
                out = conjugate(out, kron(z, id))

            pout = np.trace(out).real
        
            pid = np.trace(out * bell).real
            perror = pout - pid
            pxerror = np.trace(out * conjugate(bell, kron(id, x))).real
            pyerror = np.trace(out * conjugate(bell, kron(id, y))).real
            pzerror = np.trace(out * conjugate(bell, kron(id, z))).real

            totalpout += pout
            totalperror += perror
            totalpxerror += pxerror
            totalpyerror += pyerror
            totalpzerror += pzerror

    if totalpout <= 1e-12:
        totalpout = 0
    if totalperror <= 1e-12:
        totalperror = 0
    if totalpxerror <= 1e-12:
        totalpxerror = 0
    if totalpyerror <= 1e-12:
        totalpyerror = 0
    if totalpzerror <= 1e-12:
        totalpzerror = 0
    return totalpout, totalperror, totalpxerror, totalpyerror, totalpzerror
            
    

preject = [[[0 for _ in range(6)] for _ in range(6)] for _ in range(6)]
prejectx = 0
prejecty = 0
prejectz = 0
#probability of different types of errors given in terms of polynomial in P_X, P_Y, P_Z
xs = [[[0 for _ in range(6)] for _ in range(6)] for _ in range(6)]
ys = [[[0 for _ in range(6)] for _ in range(6)] for _ in range(6)]
zs = [[[0 for _ in range(6)] for _ in range(6)] for _ in range(6)]
total_count = [0, 0, 0, 0, 0, 0]

# Modeling of input errors. 5 qubits, 3 types of errors (+ no error) = 4^5 = 1024 input errors
for error in itertools.product(["I", "X", "Y", "Z"], repeat=5):
    input_error_name = "".join(error)
    pout, perror, pxerror, pyerror, pzerror  = stabiliser_error_simulation(5, 1, ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], error, fivequbit_decoder)
    # Number of errors that happen on the input
    input_error_size = 5 - input_error_name.count("I")
    preject[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += 1-pout
    total_count[input_error_size - 1] +=  perror
    xs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pxerror
    ys[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pyerror
    zs[input_error_name.count("X")][input_error_name.count("Y")][input_error_name.count("Z")] += pzerror

print(preject)
print(total_count)

for xerr in range(5):
    for yerr in range(5):
        for zerr in range(5):
            if preject[xerr][yerr][zerr] > 1e-8:
                print(f"Rejections given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {preject[xerr][yerr][zerr]}")

for xerr in range(5):
    for yerr in range(5):
        for zerr in range(5):
            if xs[xerr][yerr][zerr] > 1e-8:
                print(f"X error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {xs[xerr][yerr][zerr]}")

for xerr in range(5):
    for yerr in range(5):
        for zerr in range(5):
            if ys[xerr][yerr][zerr] > 1e-8:
                print(f"Y error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {ys[xerr][yerr][zerr]}")

for xerr in range(5):
    for yerr in range(5):
        for zerr in range(5):
            if zs[xerr][yerr][zerr] > 1e-8:
                print(f"Z error given number of input errors X: {xerr}, Y: {yerr}, Z: {zerr}, Count: {zs[xerr][yerr][zerr]}")