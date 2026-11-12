
class perfect_fast:
    duration = 3
    qubits = 17
    outputs = 1
    def preject(self, p, px, pz, py):
        return 5 * pz + 5 * py + 5 * px
    def pxerror(self, p, px, pz, py):
        return 5 * px * pz * pz + 5 * px * py * py
    def pyerror(self, p, px, pz, py):
        return 5 * py * pz * pz + 5 * py * px * px 
    def pzerror(self, p, px, pz, py):
        return 5 * pz * py * py + 5 * pz * px * px

class iceberg_fast_4:
    duration = 2
    qubits = 19
    outputs = 2
    def preject(self, p, px, pz, py):
        return  4 * px + 4 * py + 4 * pz
    def pxerror(self, p, px, pz, py):
        return  4 * py * py + 8 * px * px
    def pyerror(self, p, px, pz, py):
        return 4 * py * py
    def pzerror(self, p, px, pz, py):
        return 4 * py * py + 8 * pz * pz

class iceberg_fast_6:
    duration = 2
    qubits = 31
    outputs = 4
    def preject(self, p, px, pz, py):
        return 6 * pz + 6 * py + 6 * px
    def pxerror(self, p, px, pz, py):
        return 16 * py * py + 32 * px * px
    def pyerror(self, p, px, pz, py):
        return 16 * py * py
    def pzerror(self, p, px, pz, py):
        return 16 * py * py + 32 * pz * pz

class iceberg_fast_8:
    duration = 2
    qubits = 43
    outputs = 6
    def preject(self, p, px, pz, py):
        return 8 * pz + 8 * py + 8 * px
    def pxerror(self, p, px, pz, py):
        return 
    def pyerror(self, p, px, pz, py):
        return 
    def pzerror(self, p, px, pz, py):
        return 

class zrep:
    duration = 2
    qubits = 2
    outputs = 1
    def preject(self, p, px, pz, py):
        return 2 * py + 2 * px
    def pxerror(self, p, px, pz, py):
        return py * py + px * px
    def pyerror(self, p, px, pz, py):
        return 2 * px * py
    def pzerror(self, p, px, pz, py):
        return 2 * pz

class xrep:
    duration = 2
    qubits = 2
    outputs = 1
    def preject(self, p, px, pz, py):
        return 2 * py + 2 * pz
    def pxerror(self, p, px, pz, py):
        return 2 * px
    def pyerror(self, p, px, pz, py):
        return 2 * pz * py
    def pzerror(self, p, px, pz, py):
        return pz * pz + py * py