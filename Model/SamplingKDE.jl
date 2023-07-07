module SamplingKDE
    """
    Module to sample from the soil moisture-evap KDE.
    """
    
    using PyCall
    using Distributions
    export sampling

    function sampling(pds_per_evap_soil, evap)
        py"""
        import numpy as np
        if evap == True:
            vector = np.load("/Model/KDE/evap_range_3D_complete.npy")
        else:
            vector = np.load("/Model/KDE/runoff_range_3D_complete.npy")
        def sample(pds_per_evap_soil):
            samples = np.random.choice(vector, p=pds_per_evap_soil)
            return samples
        """
        samples = py"sample"(pds_per_evap_soil, evap)
    end

    function sampling_jl(pds_per_evap_soil, evap_vector)
        """
        Function to sample evaporation from the probability density (for fixed value of the soil moisture and ssr).
        """

        sampled = 0
        if isprobvec(pds_per_evap_soil) == false
            print(pds_per_evap_soil, sum(pds_per_evap_soil))
        end
        a = Categorical(pds_per_evap_soil)
        DiscreteNonParametric{Int64,Float64,Base.OneTo{Int64},Array{Float64,1}}(Base.OneTo(50), pds_per_evap_soil)
        for i = 1:100
            sampled += evap_vector[rand(a)]
        end
        return sampled/100

    end
end
