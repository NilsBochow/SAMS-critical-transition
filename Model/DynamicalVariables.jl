module DynamicalVariables
    """
    Module with definition of variables: P, E, W, A, R and Q.
    """
    
    include("../Seasonality-ERA5/Constants.jl")
    include("../Seasonality-ERA5/SamplingKDE.jl")
    import .Constants; const C = Constants
    using .SamplingKDE
    using Statistics

    export fillP, fillR, Pi2, fillQ, fillW2, ma_t, ms_t, fillE_sampled

    function P_box(ma_value, p0, p1, p2, ath0, ath1, i, t)
        """
        Calculate the precipitation P depending on the box i.
        """
        p0 = p0[i]
        p1 = p1[i]
        p2 = p2[i]
        P = p0*exp(p1*ma_value)+p2
        if P<0
            P=0
        end
        return P
    end

    function fillP(ma_vector, p0, p1, p2, n, ath0, ath1, t)
        """
        Calculate P (precipitation) in every box.
        """
        P_vector = zero(ma_vector)
        for i = 1 : n
            P_vector[i]=P_box(ma_vector[i], p0, p1, p2, ath0, ath1, i, t)
        end
        return P_vector
    end

    function E_stat_ssr_soil(KDE_vector, evap_vector)
        """
        Calling the function to sample E from the KDE.
        """
        E = SamplingKDE.sampling_jl(KDE_vector, evap_vector)
        return E
    end

    function fillE_sampled(KDE_vector, evap_vector, n, tlen)
        """
        Calculate E (evaporation) in every box.
        """
        E_vector = fill(0.0, n)
        for i= 1:n
            E_vector[i] = E_stat_ssr_soil(KDE_vector[i, :], evap_vector)
        end
        return E_vector
    end

    function R_box(ms_value, r1, r2, r3, sth, i)
        """
        Calculate R (runoff).
        """
        r1 = r1[i]
        r2 = r2[i]
        r3 = r3[i]
        if ms_value > 0 && ms_value < sth
            R = exp(r1 * ms_value - r2) + r3
        elseif ms_value > sth
            R = ms_value - sth + exp(r1 * ms_value -r2) + r3
        elseif ms_value <= 0
            R = 0
        else
            println(ms_value)
        end
        return R
    end

    function fillR(ms_vector, r1, r2, r3, sth, n, tlen)
        """
        Calculate R (runoff) for every box.
        """
        R_vector = zero(ms_vector)
        for i = 1 : n
            R_vector[i]=R_box(ms_vector[i], r1, r2, r3, sth, i)
        end
        return R_vector
    end

    function Q(P_value, SH_value, RH_value)
        Q = P_value * 2600000.0 / 3600.0 + SH_value + RH_value
        return Q
    end

    function Pi2(Q_mat, pi, t, Q_o)
        pi = mean(Q_mat[11 : 91, t]) - Q_o
        if pi<0
            pi = 0
        end
        return pi
    end

    function fillQ(P_vector, SH_vector, RH_vector, n, tlen)
        Q_vector = zero(P_vector)
        for i= 1 : n
            Q_vector[i]=Q(P_vector[i], SH_vector[i], RH_vector)
        end
        return Q_vector
    end

    function W(i,w0,L,pi,W_vec, t)
        """
        Calculate the winds.
        """

        W = (16.5 - w0) * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))  +
            w0 * L * pi * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4))) #+(i-1)^2/3000
        if L*pi>1
            W = (16.5 - w0) * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))  +
                w0 *1* (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))#+(i-1)^2/3000

        end

        if i > 70
            W = (16.5 - w0) * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))  +
                w0 * L * pi * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4))) +(i-70)^2/400*wind_cos_value_quadratic(t)
            if L*pi>1
                W = (16.5 - w0) * (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))  +
                    w0 *1* (1.0 + 1.0 / (1.0 +  exp((.06 + wind_cos_value(t)*1)* (i-1) - 3.4)))+(i-70)^2/400*wind_cos_value_quadratic(t)

            end
        end



        return W
    end

    function wind_cos_value(t)
        1*(cos(t*7.17186433e-04-3.1415926535897)*-0.03 - 0.03)
    end

    function wind_cos_value_quadratic(t)
        1*(cos(t*7.17186433e-04-3.1415926535897)*0.5 + 0.5)
    end

    function fillW2(n, w0, L, pi, t)
        W_vec = Array{Float64}(undef, n)
        for i = 1 : n
            W_vec[i] = W(i, w0, L, pi, W_vec, t)
        end
        replace!(x -> x>30 ? 30 : x, W_vec)
        replace!(x -> x<2 ? 2 : x, W_vec)
        return W_vec
    end

    function ma_t(t, n, E_mat, P_mat, W_mat, ma, l, noise)
        """
        Integrate the atmospheric moisture ma.
        """
        for i = 2:n
            ma[i, t + 1] = ma[i, t] + E_mat[i, t] - P_mat[i, t] +
                (W_mat[i - 1, t]* ma[i - 1, t] - W_mat[i, t] * ma[i, t]) / (l) + noise[i, t]
            if ma[i, t + 1] < 0
                ma[i, t + 1] = 0
            end
        end
        return ma[:, t + 1]
    end

    function ms_t(t, n, E_mat, P_mat, R_mat, ms)
        """
        Integrate the soil moisture ms.
        """
        for i=2:n
            ms[i, t + 1] = ms[i, t] +  P_mat[i,t] - E_mat[i,t]  - R_mat[i,t]
            if ms[i, t + 1] > 1500
                ms[i, t + 1] = 1500
            end
            if ms[i, t + 1] < 35
                ms[i, t + 1] = 35
            end
        end
        return ms[:, t + 1]
    end

end
