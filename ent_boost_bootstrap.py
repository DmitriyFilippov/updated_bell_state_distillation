
from util import NoiseConfiguration
from entanglement_boosting import (perform_distillation, CircuitWithAdditionalProperties, perform_parallel_simulation)
from Factory import (Factory, Hardware, GateSurfaceCode)
import csv

num_shots = 1000000000
parallelism = 8

def run(hardware_error, bell_error):

    hw = Hardware(hardware_error, bell_error, 1e-6, 50e-9, GateSurfaceCode) # type: ignore
    

    noise_conf = NoiseConfiguration(
            single_qubit_gate_error_probability=hw.phys_e,
            two_qubit_gate_error_probability=hw.phys_e,
            reset_error_probability=hw.phys_e,
            measurement_error_probability=hw.phys_e,
            idle_error_probability=hw.phys_e)

    with open(f"entanglement_boosting-{hw.phys_e}-{hw.bell_e}-v2.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        for bell_distance in range(3,6,2):
            for surface_distance in range(bell_distance, 8,2):
                for discard_rate in [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
                    circuits: CircuitWithAdditionalProperties = perform_distillation(
                        bell_distance, surface_distance, noise_conf, hw.bell_e, True)
                    results = perform_parallel_simulation(
                        circuits, num_shots, parallelism, 2 ** 20, False)
                    num_unconditionally_discarded = results.num_discarded_samples
                    num_total_samples = len(results)
                    num_valid = 0
                    num_wrong = 0
                    num_discarded = num_unconditionally_discarded
                    discarding = True
                    gap_threshold = results.max_gap()
                    for (gap, bucket) in enumerate(results.buckets):
                        if discarding and len(bucket) + num_discarded <= num_total_samples * discard_rate:
                            num_discarded += len(bucket)
                            continue

                        if discarding:
                            discarding = False
                            gap_threshold = gap

                        num_valid += bucket.num_valid_samples
                        num_wrong += bucket.num_wrong_samples

                    factory = Factory(name = f"Entanglement Boosting {bell_distance} {surface_distance} {discard_rate}", 
                                    pout = (num_valid + num_wrong) / num_shots,
                                    pxerror = (num_wrong / (num_valid + num_wrong)) / 3,
                                    pyerror = (num_wrong / (num_valid + num_wrong)) / 3,
                                    pzerror = (num_wrong / (num_valid + num_wrong)) / 3,
                                    duration = 1,
                                    qubits = 2 * (surface_distance**2) - 1,
                                    inputs = bell_distance ** 2,
                                    outputs = 1)
                    writer.writerow([factory.name, factory.pout,factory.pxerror,factory.pzerror,factory.pyerror, factory.duration,factory.qubits, factory.inputs, factory.outputs])


def main() -> None:
    run(1e-4, 1e-2)
    run(1e-3, 1e-2)
    run(1e-4, 5e-2)
    run(1e-3, 5e-2)

if __name__ == '__main__':
    main()
