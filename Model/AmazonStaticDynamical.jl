module AmazonStaticDynamical
    include("../Seasonality-ERA5/DynamicalVariables.jl")
    include("../Seasonality-ERA5/Constants.jl")
    #include("../Seasonality-ERA5/InitialMoisture.jl")
    include("../Seasonality-ERA5/MoistureDrought.jl")
    #include("../Seasonality-ERA5/InitialSSR.jl")
    include("../Seasonality-ERA5/LoadEvapKDE.jl")
    using .DynamicalVariables
    using .LoadEvapKDE
    #using .InitialMoisture
    #using .InitialSSR
    using .MoistureDrought
    using JLD
    import .Constants; const C = Constants
    export amazon_static, amazon_dynamical


    function return_undef_arrays()
        """
        Initialize arrays for saving.
        """

        W_mat = initialize_array()
        R_mat = similar(W_mat)
        P_mat = similar(W_mat)
        E_mat = similar(W_mat)
        Q_mat = similar(W_mat)
        ms = similar(W_mat)
        ma = similar(W_mat)
        ssr = Array{Float64,2}(undef, 11, C.tlen_dyn)
        noise = C.std * randn(C.n, C.tlen_dyn)
        SH = ones(C.n, C.tlen_dyn)
        pi = 0
        return W_mat, R_mat, P_mat, E_mat, Q_mat, ms, ma, ssr, noise, SH, pi
    end

    function initialize_array()
        """
        Function to initialize arrays with dimension n x tlen
        """
        
        array_n_tlen = Array{Float64,2}(undef, C.n, C.tlen_dyn)
        return array_n_tlen
    end

    function initialize_moisture(ma, ms, vapor_init, soil_init)
        """
        Initial values for atmospheric moisture and soil moisture.
        Atm. moisture initial values in Box 1 taken from ERA5-reanalyis.
        """
        ma[:, 1] .= C.ma0
        ma[1, :] .= vapor_init[1:C.tlen_dyn]
        ms[:, 1] .= 921
        ms[1, :] .= 921#soil_init[1:C.tlen_dyn]
        return ma, ms
    end

    function initialize_ssr(ssr_init)
        """
        Function to select the surface solar radiation values over Amazonia/trajectory for all time steps of the model.
        """
        ssr = ssr_init[:, 1:C.tlen_dyn]
    end
    function intialize_RH()
        """
        Function to approximate radiative heating over Amazonia/trajectory.
        """
        RH = cos.(range(0,length= C.tlen_dyn, stop =floor(C.tlen_dyn/(24*365/12))) .* (0.52354912).+0.96716066)*24.6 .-90.1
    end
    function initialize_KDE_arrays(KDE_3D)
        """
        Initializing arrays for the sampling from the KDE.
        """
        soil_value = zeros(Int16, C.n)
        ssr_value = zeros(Int16, C.n)
        KDE_2D = zeros(C.n, size(KDE_3D)[2])
        return soil_value, ssr_value, KDE_2D
    end

    function get_KDE_2D(soil_range, ssr_range, ms, ssr, KDE_3D, t)
        """
        Function to get probabilty densities for sepcific ssr and soil moisture from the 3D KDE at every box for time t.
        """
        soil_value, ssr_value, KDE_2D = initialize_KDE_arrays(KDE_3D)
        for i = 1:C.n
            # Find index in soil and ssr vectors where soil/ssr value are equal for every box
            soil_value[i] = convert(Int16, searchsortedfirst(soil_range, ms[i, t]) - 1)
            ssr_value[i] = convert(Int16, searchsortedfirst(ssr_range, ssr[i, t]) - 1)
            if soil_value[i]==0
                println(ms[i,t])
            end
            # select values from 3D KDE
            KDE_2D[i, :] = KDE_3D[ssr_value[i], :, soil_value[i]]
        end
        return soil_value, ssr_value, KDE_2D
    end

    function initialize_deforestation_arrays()
        e_REW = zeros(C.n)
        e_integrated_drought = ones(C.n)
        soil_baseline = similar(e_REW)

        return e_REW, e_integrated_drought, soil_baseline
    end

    function manual_deforestation!(stability, SH, defo0, defo1, sh1)
        """
        Function to deforest the trajectory. Changes the evapotranspiration and sensible heating when box gets deforested.
        """
        println(defo1)
        defo_factor = ones(C.n, C.tlen_dyn)
        if stability == true
            defo_factor[1: defo0, :] .= 0.6
            SH[1: defo0, :] .= C.sh0 * (1. + sh1)
        else
            for i in defo0 + 1 : defo1
                defo_factor[i, convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
                    C.tlen_dyn / 2)) + 1 : end] .= 0.6
                SH[i, convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
                    C.tlen_dyn / 2)) + 1 : end] .= C.sh0 * (1. + sh1)
            end
        end
        return defo_factor, SH
    end

    function manual_reforestation2!(SH, defo_factor, defo0, defo1, sh1)
        """
        Function for reforestation experiments. Basically, inversion of the manual_deforestation function.
        """
        defo_factor .= 0.6
        SH .= C.sh0 * (1. + sh1)
        for i in defo0 + 1 : defo1
            #defo_factor[i, 1 : convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
            #    C.tlen_dyn / 2)) + 1 ] .= 0.6
            defo_factor[(102-i), convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
                C.tlen_dyn / 2)) + 1 : end] .=1
            #SH[i, 1 : convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
            #    C.tlen_dyn / 2)) + 1 ] .= C.sh0 * (1. + sh1)
            SH[(102-i), convert(Int64, floor(C.tlen_dyn / 4 + (i-1) / (C.n-1) *
                C.tlen_dyn / 2)) + 1 : end] .= C.sh0
        end
        return defo_factor, SH
    end

    function soil_baseline_save(t, ms, soil_baseline)
        """
        Saving the soil moisture baseline from model run without deforestation to calculate the soil moisture deficit
        when reforesting the rainforest.
        """
        soil_baseline[:] = minimum(ms[:, 4000:convert(UInt32, floor(C.tlen_dyn/4))], dims = 2)
        #save("./Model/saved_variables/soil_baseline.jld", "soil_baseline", soil_baseline)
        return soil_baseline
    end
    function drought_effects(t, ms, defo_factor, e_REW, e_integrated_drought, soil_baseline, integrated_soil_deficit, SH_drought, sh1)
        """
        Functions to call the functions for the soil-moisture-vegetation feedback.
        Checks if box is already deforested in every time step. If it is deforested due to the SMD, the sensible heat increases.
        If the box is deforested due to manual deforestation but the SMD does not exceed the threshold, the box is deforested nonetheless.
        """

        if t > C.tlen_dyn/4+1
            integrated_soil_deficit = integrated_drought(integrated_soil_deficit, ms[:, t], soil_baseline)
            e_integrated_drought = check_integrated_drought(integrated_soil_deficit, 0.6)
        else
            e_integrated_drought = ones(C.n)
        end
        for i=1:C.n
            if e_integrated_drought[i] == 0.6
                SH_drought[i] = C.sh0 * (1 + sh1)
            end
            if defo_factor[i, t] == 0.6 && e_integrated_drought[i] != 0.6
                e_integrated_drought[i] = defo_factor[i, t]
                SH_drought[i] = C.sh0 * (1 + sh1)
            end
        end
        return e_integrated_drought, SH_drought
    end

    function amazon_static(L, w0, defo0, defo1, dr, Q_o, vapor_init, soil_init, ssr_init, KDE_3D, ssr_range, evap_vector, soil_range)
        """
        Function that calls all the functions for calculations without deforestation, used for dimensionality factor L.
        """
        W_mat, R_mat, P_mat, E_mat, Q_mat, ms, ma, ssr, noise, SH, pi  = return_undef_arrays()
        SH = SH * C.sh0
        RH = intialize_RH()
        ma, ms = initialize_moisture(ma, ms, vapor_init, soil_init)
        ssr = initialize_ssr(ssr_init)

            @inbounds for t = 1 : C.tlen_dyn-1
                #println(t)
                P_mat[:, t] = fillP(ma[:, t], C.p0, C.p1, C.p2, C.n, C.ath0, C.ath1, t)
                soil_value, ssr_value, KDE_2D = get_KDE_2D(soil_range, ssr_range, ms, ssr, KDE_3D, t)

                E_mat[:, t] = fillE_sampled(KDE_2D, evap_vector, C.n, C.tlen_dyn)
                R_mat[:, t] = fillR(ms[:, t], C.r1, C.r2, C.r3, C.sth, C.n, C.tlen_dyn)
                Q_mat[:, t] = fillQ(P_mat[:, t], SH[:, t], RH[t], C.n, C.tlen_dyn)
                pi = Pi2(Q_mat, pi, t, Q_o[t])
                W_mat[:, t] = fillW2(C.n, w0, L, pi, t)
                ma[:, t + 1] = ma_t(t, C.n, E_mat, P_mat, W_mat, ma, C.l, noise)
                ms[:, t + 1] = ms_t(t, C.n, E_mat, P_mat, R_mat, ms)

            end
        return ma[:, :], Q_mat[:, :], P_mat[:, :], W_mat[:, :], E_mat[:, :], ms[:, :], R_mat[:, :]
    end

    function amazon_dynamical(L, w0, defo0, defo1, Q_o, refo, sh1, vapor_init, soil_init, ssr_init, KDE_3D, ssr_range, evap_vector, soil_range, stability, soil_baseline, integrated_soil_deficit)
        """
        Function that calls all the functions for calculations.
        """
        W_mat, R_mat, P_mat, E_mat, Q_mat, ms, ma, ssr, noise, SH, pi  = return_undef_arrays()
        SH = SH * C.sh0
        RH = intialize_RH()
        ma, ms = initialize_moisture(ma, ms, vapor_init, soil_init)
        ssr = initialize_ssr(ssr_init)
        e_REW, e_integrated_drought, soil_baseline = initialize_deforestation_arrays()
        defo_factor, SH = manual_deforestation!(stability, SH, defo0, defo1, sh1)
        if refo == 1
            defo_factor, SH = manual_reforestation2!(SH, defo_factor, defo0, defo1, sh1)
            RH = RH[end:-1:1]
            Q_o = Q_o[:, end:-1:1]
            ssr = ssr[:, end:-1:1]
            ma[1, :] = ma[1, end:-1:1]
            ms[1, :] = ms[1, end:-1:1]
        end

        @inbounds for t = 1 : C.tlen_dyn - 1

            P_mat[:, t] = fillP(ma[:, t], C.p0, C.p1, C.p2, C.n, C.ath0, C.ath1, t)
            soil_value, ssr_value, KDE_2D = get_KDE_2D(soil_range, ssr_range, ms, ssr, KDE_3D, t)
            if t == convert(UInt32, floor(C.tlen_dyn/4))
                if refo != 1
                    soil_baseline = soil_baseline_save(t, ms, soil_baseline)
                else
                    #soil_baseline = load(pwd()*"./Model/saved_variables/soil_baseline.jld", "soil_baseline")
                    soil_baseline = soil_baseline
                end
            end
            e_integrated_drought, SH[:, t] = drought_effects(t, ms, defo_factor, e_REW, e_integrated_drought, soil_baseline, integrated_soil_deficit, SH[:, t], sh1)
            E_mat[:, t] = fillE_sampled(KDE_2D, evap_vector, C.n, C.tlen_dyn).*e_integrated_drought
            R_mat[:, t] = fillR(ms[:, t], C.r1, C.r2, C.r3, C.sth, C.n, C.tlen_dyn)
            Q_mat[:, t] = fillQ(P_mat[:, t], SH[:, t], RH[t], C.n, C.tlen_dyn)
            pi = Pi2(Q_mat, pi, t, Q_o[t])
            W_mat[:, t] = fillW2(C.n, w0, L, pi, t)
            ma[:, t + 1] = ma_t(t, C.n, E_mat, P_mat, W_mat, ma, C.l, noise)
            ms[:, t + 1] = ms_t(t, C.n, E_mat, P_mat, R_mat, ms)
        end
        return ma[:, :], Q_mat[:, :], P_mat[:, :], W_mat[:, :], E_mat[:, :], ms[:, :], R_mat[:, :], soil_baseline[:]
    end

end
