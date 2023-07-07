module Constants
    """
    Module to define all the constants/parameters used in the model.
    """

    include("../Seasonality-ERA5/InterpolateParameters.jl")
    using .InterpolateParameters
    const n = 101
    const year = 40
    const tlen_stat= 365*year*24
    const tlen_dyn= 365*year*24

    const fb = 96
    const ath0=35
    const ath1=55
    const std = 0.
    const ma0 = 53.8703304702
    const ma1 = 52
    const ms0 = 120
    const ms1 = 120
    const rng = 84476
    const l = 30

    
    const e0 = 0.1614321434117725
    const e1 = -0.9354672249495436
    const e2 = 13.818289560334762



    const sth1=110
    const sth=1450


    
    r1_array = ones(11)*0.02356
    for i in [1, 5, 6, 7, 8, 9, 11]
        r1_array[i] =0.01456
    end
    r1_array[4]=0.02156
    r1_array[5]=0.01356
    r1_array[6]=0.01696
    r1_array[9]=0.01396
    r1_array[end]=0.01316
    r2_array = ones(11)*33.794078433062886
    for i in  [1, 5, 6, 7, 8, 9, 11]
        r2_array[i] = 19.48
    end
    r2_array[6] = 22.48
    r3_array = [ 0.004, 0.004, 0.004,
        0.004, 0.004, 0.004, 0.004,
        0.004, 0.004, 0.004, 0.004]
    r3_array = ones(11)*0.0018683905259895554
    for i in  [1, 5, 6, 7, 8, 9, 11]
        r3_array[i] = 0.0013383905259895554
    end


    r1 = InterpolateParameters.interp(r1_array)
    r2 = InterpolateParameters.interp(r2_array)
    r3 = InterpolateParameters.interp(r3_array)

    p0_array = ones(11)*0.07920957028545568
    p1_array = ones(11)*0.0426946155306704
    p2_array = ones(11)*-0.40701849323166366

    p0 = InterpolateParameters.interp(p0_array)
    p1 = InterpolateParameters.interp(p1_array)
    p2 = InterpolateParameters.interp(p2_array)

    const s0 = 0.028871004086204754
    const s1 =  0.0007981567618544439
    const sh0 = 24.84022197313726
    const sh1 = 0.44
    const rh0 = -97.4769925204
    const emax = 0.16

end
