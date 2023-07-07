module InterpolateParameters
    """
    Module to interpolate values at every model box (ERA5 trajectory only consists of 10 boxes, 1 ERA5 box = 10 model boxes).
    """
    
    using DifferentialEquations
    using Random
    using PyCall

    export readNC


    function interp(array)
        py"""
        import numpy as np
        from matplotlib import pyplot as plt

        def interpolate(array):
            monthlyArray = np.arange(0, 110, 10)

            array_interp = np.interp(np.arange(101), monthlyArray, array)
            return array_interp
        """
        x = py"interpolate"(array)
    end
end
