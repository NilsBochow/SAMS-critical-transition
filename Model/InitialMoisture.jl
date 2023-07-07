module InitialMoisture
    """
    Module to convert, interpolate (at every hour) and save initial atmospheric moisture (vapor) and soil moisture from ERA5-reanalysis from .npy to .jld.
    """

    using DifferentialEquations
    using Random
    using PyCall
    using JLD

    export initMoisture, generatecos, initSoil, initQ_o

    function initMoisture(year)
        py"""
        import numpy as np
        from matplotlib import pyplot as plt
        import os

        print(os.getcwd())
        def realData(year):
            vapor = np.load(os.getcwd() + '/saved_variables/vapor_SA.npy')
            indices = np.load(os.getcwd() + "/saved_variables/indices_ERA5.npy", allow_pickle = True)
            vapor = np.mean(vapor.reshape((492, 221 * 281))[:12*year][:, indices[-1]], axis = 1).flatten()
            return vapor

        def interpolate(year):
            vapor = realData(year)
            monthlyArray = np.arange(0, year*24*365, 365/12*24)
            vapor_interp = np.interp(np.arange(year*24*365), monthlyArray, vapor)
            return vapor_interp
        """
        x = py"interpolate"(year)
        save(pwd() * "/saved_variables/vapor_init.jld", "vapor_init", x)
        return x
    end

    function initSoil(year)
        py"""
        import numpy as np
        from matplotlib import pyplot as plt
        import os

        def realData(year):
            soil = np.load(os.getcwd() + '/Model/saved_variables/soil_deep_SA.npy')
            indices = np.load(os.getcwd() + "/Model/saved_variables/indices_ERA5.npy", allow_pickle = True)
            soil = np.mean(soil.reshape((492, 221 * 281))[:12*year][:, indices[-1]], axis = 1).flatten()
            return soil

        def interpolate(year):
            soil = realData(year)
            monthlyArray = np.arange(0, year*24*365, 365/12*24)
            soil_interp = np.interp(np.arange(year*24*365), monthlyArray, soil)
            return soil_interp
        """
        x = py"interpolate"(year)
        #println(x)
        save(pwd() * "/Model/saved_variables/soil_init.jld", "soil_init", x)
        return x
    end

    function initQ_o(year)
        py"""
        import numpy as np
        from matplotlib import pyplot as plt

        def realData(year):
            Q_o = np.load(os.getcwd() + '/Model/saved_variables/tot_heat_Q_o.npy')
            Q_o = Q_o[:12*year].flatten()
            return Q_o

        def interpolate(year):
            Q_o = realData(year)
            monthlyArray = np.arange(0, year*24*365, 365/12*24)
            Q_o_interp = np.interp(np.arange(year*24*365), monthlyArray, Q_o)
            return Q_o_interp
        """
        x = py"interpolate"(year)
    end

    function generatecos(tlen)
        Random.seed!(123)
        W = rand!(zeros(10))/10
        V = rand(0:0.01:3, 10)
        cos2 = (64/30*cos.(range(0,length= tlen, stop =floor(tlen/(24*365/12))-1) .* (2*3.14159/12) .-1)).+48
        for i=1:5
          cos2 += (V[i]*cos.(range(0,length= tlen, stop =floor(tlen/(24*365/12))-1) .* (2*3.14159/12 .+W[i]) .-1))
        end
        for i = 5:10
          cos2 += (V[i]/2*sin.(range(0,length= tlen, stop =floor(tlen/(24*365/12))-1) .* (2*3.14159/12 .+10*W[i]) .-1))
        end
        return cos2
    end

end
