module CalculateDynamical
    """
    Module/File to execute, that calls all relevant functions and saves the solution arrays.
    Every variable at every time step (hourly) in every box for all simulation parameters
    gets saved, so the solution arrays are large.
    Change the corresponding saving functions if you do not want the variables at every time step.
    """

    println(Threads.nthreads())
    include("../Seasonality-ERA5/AmazonStaticDynamical.jl")
    include("../Seasonality-ERA5/PlotStatic.jl")
    include("../Seasonality-ERA5/Constants.jl")
    include("../Seasonality-ERA5/CalculateStatic.jl")
    include("../Seasonality-ERA5/PlotDynamical.jl")
    include("../Seasonality-ERA5/InitialMoisture.jl")
    include("../Seasonality-ERA5/LoadEvapKDE.jl")
    using .LoadEvapKDE
    using .AmazonStaticDynamical
    import .Constants; const C = Constants
    import .CalculateStatic
    import .PlotDynamical
    using .InitialMoisture

    using PyCall
    using Plots
    import PyPlot; const plt = PyPlot
    pygui(true)
    using LaTeXStrings
    using Statistics
    using JLD

    function initalize_Q_o_array(Q_o_amp, Q_o)
        """
        Function that approximates the heating over the Atlantic ocean as cosine function.
        """
        Q_o_array = cos.(range(0,length= C.tlen_dyn, stop =floor(C.tlen_dyn/(24*365/12))) .* (0.52354912).+1.02160148)*-Q_o_amp .+Q_o
    end

    function initialize_arrays_static(Q_o_amp_Range, Q_oRange)
        """
        Function to initialize arrays for the static run.
        """
        ma = zeros(C.n, C.tlen_dyn, size(Q_o_amp_Range)[1], size(Q_oRange)[1])
        ms = similar(ma)
        Q_mat = similar(ma)
        P_mat = similar(ma)
        R_mat = similar(ma)
        W_mat = similar(ma)
        E_mat = similar(ma)
        return ma, ms, Q_mat, P_mat, R_mat, W_mat, E_mat
    end

    function initialize_arrays_dynamic(Q_o_amp_Range, Q_oRange, w0Range)
        """
        Function to initialize arrays for the dynamical run.
        """
        ma = zeros(C.n, C.tlen_dyn, size(Q_o_amp_Range)[1], size(Q_oRange)[1], size(w0Range)[1])
        soil_baseline = zeros(C.n, size(Q_o_amp_Range)[1], size(Q_oRange)[1], size(w0Range)[1])
        ms = similar(ma)
        Q_mat = similar(ma)
        P_mat = similar(ma)
        R_mat = similar(ma)
        W_mat = similar(ma)
        E_mat = similar(ma)
        return ma, ms, Q_mat, P_mat, R_mat, W_mat, E_mat, soil_baseline
    end

    function set_sh1(defo1)
        """
        Function to set the sensible heating
        """
        if defo1 == 0
            sh1 = 0
        else
            sh1 = 0.44
        end
        return sh1
    end


    function save_variables(defo_string, ma, P_mat, W_mat, E_mat, R_mat, ms, Q_mat, soil_baseline_saved)
        """
        Function to save all variable arrays as .jld (hdf5 format).
        """
        save(pwd() * "/saved_variables/"*defo_string*"/moisture.jld", "ma", ma)
        save(pwd() * "/saved_variables/"*defo_string*"/precip.jld", "P_mat", P_mat)
        save(pwd() * "/saved_variables/"*defo_string*"/wind.jld", "W_mat", W_mat)
        save(pwd() * "/saved_variables/"*defo_string*"/evap.jld", "E_mat", E_mat)
        save(pwd() * "/saved_variables/"*defo_string*"/runoff.jld", "R_mat", R_mat)
        save(pwd() * "/saved_variables/"*defo_string*"/soil.jld", "ms", ms)
        save(pwd() * "/saved_variables/"*defo_string*"/heat.jld", "Q_mat", Q_mat)
        save(pwd() * "/saved_variables/"*defo_string*"/soil_baseline.jld", "soil_baseline", soil_baseline_saved)
    end

    function callingDynamical(drRange, Q_oRange, Q_o_amp_Range, w0Range)
        """
        Function that calls the simulation for all combinations of simulation parameters. Uses Julia Multithreading.
        """
        cluster = true
        ma, ms, Q_mat, P_mat, R_mat, W_mat, E_mat, soil_baseline_saved = initialize_arrays_dynamic(Q_o_amp_Range, Q_oRange, w0Range)

        KDE_data = load(pwd()*"/saved_variables/KDE_data.jld")#
        z_norm, ssr_range, evap_vector, soil_deep = LoadEvapKDE.loadArrayKDE()

        vapor_init= load(pwd()*"/saved_variables/vapor_init.jld", "vapor_init")#
        ssr_init = transpose(load(pwd()*"/saved_variables/ssr_init.jld", "ssr"))
        soil_init = load(pwd()*"/saved_variables/soil_init.jld", "soil_init")
        println("here")
        L = calculateL(cluster, drRange, Q_oRange, Q_o_amp_Range, vapor_init, soil_init, ssr_init, z_norm, ssr_range, evap_vector, soil_deep)
        stability = false
        for refo = 0
            if refo == 1
                soil_baseline_load = load(pwd()*"/saved_variables/defo_new/soil_baseline.jld", "soil_baseline")
            else
                soil_baseline = zeros(101)
            end
            Threads.@threads for Q_o_amp in Q_o_amp_Range
                c = CalculateStatic.ThreadingIndex(Q_o_amp, Q_o_amp_Range)
                Threads.@threads for Q_o in Q_oRange
                    d = CalculateStatic.ThreadingIndex(Q_o, Q_oRange)
                    Q_o_array = initalize_Q_o_array(Q_o_amp, Q_o)
                    Threads.@threads for w0 in w0Range
                        f = CalculateStatic.ThreadingIndex(w0, w0Range)
                        integrated_soil_deficit = zeros(C.n, size(Q_o_amp_Range)[1], size(Q_oRange)[1], size(w0Range)[1])
                        if refo == 1
                            soil_baseline = soil_baseline_load[:, c, d, f]
                        else
                            nothing
                        end
                        println("Q_o: $(Q_o)", "Q_o_amp: $(Q_o_amp)", "w0: $(w0)", "Threadid: ", Threads.threadid())
                        defo1 = 101
                        sh1 = set_sh1(defo1)
                        @time begin
                            ma[:, :, c, d, f], Q_mat[:, :, c, d, f], P_mat[:, :, c, d, f], W_mat[:, :,c ,d, f], E_mat[:, :, c, d, f], ms[:,: , c, d, f], R_mat[:, :, c, d, f], soil_baseline_saved[:, c, d, f] =
                                amazon_dynamical(L[c,d], w0, 0, defo1,  Q_o_array, refo, sh1, vapor_init, soil_init, ssr_init, z_norm, ssr_range, evap_vector, soil_deep, stability, soil_baseline, integrated_soil_deficit[:, c, d, f])
                        end
                        for box = 1:10:101
                            #defo = 0
                            
                            PlotDynamical.plotTime(false, E_mat[:, :, c ,d], ms[:, :, c ,d], ma[:, :, c ,d],
                                W_mat[:, :, c ,d], P_mat[:, :, c ,d], R_mat[:, :, c ,d], Q_mat[:, :, c ,d],
                                Q_o, Q_o_amp, Q_o_array, w0, box, refo, defo1, C.tlen_dyn, "tests")
                            
                        end
                    end
                end
            end
            #save_variables("defo_drought_new", ma, P_mat, W_mat, E_mat, R_mat, ms, Q_mat, soil_baseline_saved)
            save(pwd() * "/saved_variables/"*"defo_drought_new"*"/wind5.jld", "W_mat", W_mat[40,:,:,:])
            save(pwd() * "/saved_variables/"*"defo_drought_new"*"/ma5.jld", "ma",ma[40,:,:,:])
            save(pwd() * "/saved_variables/"*"defo_drought_new"*"/precip5.jld", "P_mat", P_mat[40,:,:,:])


        end
        println("----dynamical finished----")
    end

    function calculateL(cluster, drRange, Q_oRange, Q_o_amp_Range, vapor_init, soil_init, ssr_init, z_norm, ssr_range, evap_vector, soil_deep)
        """
        Function to calculate the dimensionality factor L.
        """
        ma, ms, Q_mat, P_mat, R_mat, W_mat, E_mat = initialize_arrays_static(Q_o_amp_Range, Q_oRange)
        L = zeros(size(Q_o_amp_Range)[1], size(Q_oRange)[1])
        dr = 0.16
        Threads.@threads for Q_o_amp in Q_o_amp_Range
            c = CalculateStatic.ThreadingIndex(Q_o_amp, Q_o_amp_Range)
            Threads.@threads for Q_o in Q_oRange
                d = CalculateStatic.ThreadingIndex(Q_o, Q_oRange)
                Q_o_array = initalize_Q_o_array(Q_o_amp, Q_o)
                println("static")
                @time begin
                ma[:, :, c, d], Q_mat[:, :, c, d], P_mat[:, :, c, d], W_mat[:, :,c ,d], E_mat[:, :, c, d], ms[:,: , c, d], R_mat[:, :, c, d] = amazon_static(0, 0, 1, 0, dr, Q_o_array, vapor_init, soil_init, ssr_init, z_norm, ssr_range, evap_vector, soil_deep)
                end
                L[c, d] = 1 / maximum( mean(Q_mat[11:91,  100:end, c, d], dims = 1)[1,:] - (Q_o_array[100:end]))
                println("Lstat: ", L[c,d] )
            end
        end
        return L
    end

end



const drRange = [ .16 ]

# Simulation parameters
#const Q_oRange = [80., 100., 120., 140.]
#const Q_o_amp_Range = [20., 40.]
#const w0Range = [5.5, 8.25, 9.167, 9.9, 10.5, 11]
const Q_oRange = [120.]
const Q_o_amp_Range = [20.]
const w0Range = [9.9]


CalculateDynamical.callingDynamical(drRange, Q_oRange, Q_o_amp_Range, w0Range)
