from definitions import *
from iceberg_code_sim import iceberg_simulation



n=6
def simulator(error):
    return iceberg_simulation(error, n)


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
    pout, perror, pxerror, pyerror, pzerror, pid  = simulator(error)
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
