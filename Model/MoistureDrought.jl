module MoistureDrought
    """
    Module to calculate the integrated soil moisture deficit.
    """

    include("../Seasonality-ERA5/Constants.jl")
    include("../Seasonality-ERA5/InterpolateParameters.jl")

    using .InterpolateParameters
    import .Constants; const C = Constants
    using Interpolations

    export lengthREW, return_previous_iterable, set_e_max, integrated_drought, check_integrated_drought
    REW = zeros(C.n)
    e_max = ones(C.n)
    counter = zeros(C.n)
    iterableArray = 1:C.n
    iterableArray_previous = 1:C.n


    function return_previous_iterable()
        return iterableArray_previous
    end

    function integrated_drought(integrated_soil_deficit, soil_value, soil_baseline)
        """
        Checks if soil moisture is below baseline and if yes, the soil moisture deficit is increased.
        """

        for i = 1:C.n
            if soil_value[i] < soil_baseline[i]
                integrated_soil_deficit[i] += soil_baseline[i] - soil_value[i]
            else
                nothing
            end
        end
        return integrated_soil_deficit
    end

    function check_integrated_drought(integrated_soil_deficit, percent_e)
        """
        If integrated soil moisture deficit exceeds the set threshold, the multiplicative factor e gets set to 0.6 (corresponds to
        40% decrease of E) otherwise it stays at 1 (no decrease of E).
        """

        e = ones(C.n)
        for i = 1:C.n
            if integrated_soil_deficit[i] >450000 #if set to infinity, there is no drought induced mortality
                e[i] = percent_e
            else
                e[i] = 1
            end
        end
        return e
    end
end
