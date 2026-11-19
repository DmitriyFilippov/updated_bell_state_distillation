from Factory import *


hw = Hardware(1e-4, 1e-2, 1e-6, 50e-9, GateSurfaceCode)

def prune(factories: list[Factory], pruning_factor, min_error):
    pruned_factories: list[Factory] = []
    for factory in factories:
        flag = False
        for another_factory in factories:
            if (another_factory.spacetime * pruning_factor < factory.spacetime 
                    and ((another_factory.pxerror * pruning_factor < factory.pxerror) or (another_factory.pxerror < min_error))
                    and ((another_factory.pyerror * pruning_factor < factory.pyerror) or (another_factory.pyerror < min_error))
                    and ((another_factory.pzerror * pruning_factor < factory.pzerror) or (another_factory.pzerror < min_error))
                    and another_factory.inputs <= factory.inputs
                    and another_factory.name != factory.name): 
                flag = True 
        if not flag:
            pruned_factories.append(factory)

    return pruned_factories

