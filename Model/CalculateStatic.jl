module CalculateStatic

    include("../Seasonality-ERA5/AmazonStaticDynamical.jl")
    include("../Seasonality-ERA5/PlotStatic.jl")
    include("../Seasonality-ERA5/Constants.jl")
    using .AmazonStaticDynamical
    using .PlotStatic
    using .Constants; const C = Constants

    using Statistics
    using JLD
    using HDF5

    export callingStatic

    function ThreadingIndex(value, rangeArray)
        return findall(x->x==value, rangeArray)[1]
    end

    function callingStatic(drRange, Q_oRange, w0Range, defo1Max)
        defo0=1
        defo1=0
        L=0

        ma = zeros(C.n, C.tlen_stat, Threads.nthreads())
        Q_mat = similar(ma)
        P_mat = similar(ma)
        W_mat = similar(ma)
        E_mat = similar(ma)
        m = zeros(4, 7 , 5, 100)
        p = similar(m)
        L = calculateL(drRange, Q_oRange)
        Threads.@threads for defo1 in 0:defo1Max
            println(defo1)
            for dr in drRange
                c = ThreadingIndex(dr, drRange)
                for Q_o in Q_oRange
                    d = ThreadingIndex(Q_o, Q_oRange)
                    for w0 in w0Range
                        f = ThreadingIndex(w0, w0Range)
                        if dr == 0.16
                            sh1 = 0
                        else
                            sh1 = 0.44
                        end
                        ma[:, :, Threads.threadid()], Q_mat[:, :, Threads.threadid()], P_mat[:, :, Threads.threadid()], W_mat[:, :, Threads.threadid()], E_mat[:,:, Threads.threadid()] = amazon_static(C.n, C.tlen_stat, L[c, d], C.l, C.ma0, C.ma1, C.ms0, C.ms1, C.e0, C.e1, C.e2, C.emax, C.p0, C.p1, w0, C.ath0, C.ath1, C.sth1, C.sth, C.r1, C.r2, C.r3, defo0, defo1, dr, Q_o, C.sh0, C.sh1, C.rh0)
                        m[c, d, f, defo1+1] = ma[C.fb, C.tlen_stat-1, Threads.threadid()]
                        p[c, d, f, defo1+1] = P_mat[C.fb, C.tlen_stat-1, Threads.threadid()]
                        if defo1 in [0,50]
                            @time begin
                                println(defo1, "hier wird geplottet")
                                plotVariable!(Q_mat/10, "Q/10", C.tlen_stat, defo1)
                                plotVariable!(E_mat*100, "E", C.tlen_stat, defo1)
                                plotVariable!(ma, "ma", C.tlen_stat, defo1)
                                plotVariable!(P_mat*100, "P", C.tlen_stat, defo1)
                                savePlot("All_w0$(w0)_Qo$(Q_o)_dr$(dr)", defo1)
                                clearPlots()
                            end
                        end
                    end
                end
            end
        end
        save("Amazon_Software/Amazon Dyn/jld/amazon_stat_m.jld", "m_saved", m)
        save("Amazon_Software/Amazon Dyn/jld/amazon_stat_p.jld", "p_saved", p)
        println("----static finished----")
    end
    function calculateL(drRange, Q_oRange)
        ma = zeros(C.n, C.tlen_stat, Threads.nthreads())
        Q_mat = similar(ma)
        P_mat = similar(ma)
        W_mat = similar(ma)
        E_mat = similar(ma)
        L = zeros(size(drRange)[1], size(Q_oRange)[1])
        for dr in drRange
            c = ThreadingIndex(dr, drRange)
            for Q_o in Q_oRange
                d = ThreadingIndex(Q_o, Q_oRange)
                ma[:, :, Threads.threadid()], Q_mat[:, :, Threads.threadid()], P_mat[:, :, Threads.threadid()], W_mat[:, :, Threads.threadid()], E_mat[:,:, Threads.threadid()] = amazon_static(C.n, C.tlen_stat, 0, C.l, C.ma0, C.ma1, C.ms0, C.ms1, C.e0, C.e1, C.e2, C.emax, C.p0, C.p1, 0, C.ath0, C.ath1, C.sth1, C.sth, C.r1, C.r2, C.r3, 1, 0, dr, Q_o, C.sh0, C.sh1, C.rh0)
                L[c, d] = 1.0 / (mean(Q_mat[11 : 90, C.tlen_stat - 1, Threads.threadid()])- Q_o)
            end
        end
        return L
    end


end
