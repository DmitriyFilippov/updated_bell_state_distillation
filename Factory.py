import math
import mpmath

class ECC:
    def logical_error_rate(self, dist: int, phys_e: mpmath.mpf) -> mpmath.mpf:
        return 0
    def qubits_per_patch(self, dist: int) -> int:
        return 0
    
    def time_steps_per_error_cycle(self, dist: int) -> int:
        return 0

    def min_distance(self, operation_error_budget: mpmath.mpf, phys_e: mpmath.mpf) -> int:
        return 0
    
class GateSurfaceCode(ECC):

    def logical_error_rate(self, dist: int, phys_e: mpmath.mpf) -> mpmath.mpf:
        return 0.03 * ((phys_e / 0.01) ** int((dist+1)/2))
    
    def qubits_per_patch(self, dist: int):
        return 2 * (dist ** 2) 
    
    def time_steps_per_error_cycle(self, dist: int):
        return 6 * dist

    def min_distance(self, operation_error_budget: mpmath.mpf, phys_e: mpmath.mpf):
        min_distance = math.ceil(2 * math.log(operation_error_budget / 0.03) / math.log(phys_e / 0.01) - 1)
        if min_distance % 2 != 1:
            min_distance += 1
        return min_distance

class MeasurementSurfaceCode(ECC):

    def logical_error_rate(self, dist: int, phys_e: mpmath.mpf) -> mpmath.mpf:
        return 0.08 * ((phys_e / 0.0015) ** int((dist+1)/2))
    
    def qubits_per_patch(self, dist: int):
        return 2 * (dist ** 2)
    
    def time_steps_per_error_cycle(self, dist: int):
        return 20 * dist
    
    def min_distance(self, operation_error_budget: mpmath.mpf, phys_e: float):
        min_distance = math.ceil(2 * math.log(operation_error_budget / 0.08) / math.log(phys_e /  0.0015) - 1)
        if min_distance % 2 != 1:
            min_distance += 1
        return min_distance

class HastingsHaahCode(ECC):

    def logical_error_rate(self, dist: int, phys_e: mpmath.mpf) -> mpmath.mpf:
        return 0.07 * ((phys_e / 0.01) ** int((dist+1)/2))
    
    def qubits_per_patch(self, dist: int):
        return 4 * (dist ** 2) + 8 * (dist-1)
    
    def time_steps_per_error_cycle(self, dist: int):
        return 3 * dist
    
    def min_distance(self, operation_error_budget: mpmath.mpf, phys_e: float):
        min_distance = math.ceil(2 * math.log(operation_error_budget / 0.07) / math.log(phys_e /  0.01) - 1)
        if min_distance % 2 != 1:
            min_distance += 1
        return min_distance


class Hardware:
    def __init__(self, 
                 phys_e: mpmath.mpf, 
                 bell_e: mpmath.mpf, 
                 bell_time_in_seconds: float, 
                 op_time_in_seconds: float, 
                 ecc: ECC):
        self.phys_e = phys_e
        self.bell_e = bell_e
        self.op_time_in_seconds = op_time_in_seconds
        self.bell_time = bell_time_in_seconds / op_time_in_seconds
        self.ecc = ecc

    def get_bell_rate(self, logical_time_step):
        """
        Returns the number of bell states produced per logical time step
        """
        return logical_time_step / self.bell_time

class Factory:
    def __init__(self, name, pout: mpmath.mpf, pxerror: mpmath.mpf, pyerror: mpmath.mpf, pzerror: mpmath.mpf, duration, qubits, inputs, outputs=1):
        self.name = name 
        self.pout = pout
        self.pxerror = pxerror
        self.pyerror = pyerror
        self.pzerror = pzerror
        self.perror = pxerror + pyerror + pzerror
        self.duration = math.ceil(duration)
        self.qubits = math.ceil(qubits)
        self.spacetime = duration*qubits
        self.inputs = math.ceil(inputs)
        self.outputs = outputs
        self.production_rate = math.ceil(duration / (outputs))

