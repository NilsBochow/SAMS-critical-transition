module LoadEvapKDE
    using PyCall
    using JLD

    export loadArrayKDE, loadArrayKDE_soil_runoff

    function loadArrayKDE()
        """
        Function to load the 3D KDE and corresponding values of the SSR and evaporation obtained from ERA5 variables.
        """
        py"""
        import numpy as np
        import os
        def LoadNumpyArray():
            #KDE_ssr = np.load(os.getcwd() + "/Model/KDE/KDE_evap_ssr_soil_deep_3D.npy")
            KDE_ssr = np.load(os.getcwd() + "/KDE/KDE_evap_ssr_soil_complete_3D.npy")
            evap_vector = np.load(os.getcwd() + "/KDE/evap_range_3D.npy")
            ssr = np.load(os.getcwd() + "/KDE/ssr_range_3D.npy")
            #soil_deep = np.load(os.getcwd() + "/Model/KDE/soil_range_3D.npy")
            soil_deep = np.load(os.getcwd() + "/KDE/soil_range_3D_complete.npy")
            z_integrated = np.sum(KDE_ssr[:,:,:], axis = 1)
            z_norm = KDE_ssr[:,:,:]/z_integrated[:,None, :]
            return z_norm, ssr, evap_vector, soil_deep
        """
        z_norm, ssr, evap_vector, soil_deep = py"LoadNumpyArray"()
        save(pwd() * "/KDE/KDE_data.jld", "z_norm", z_norm, "ssr",  ssr, "evap_vector", evap_vector, "soil_deep", soil_deep)
        return z_norm, ssr, evap_vector, soil_deep
    end

    function loadArrayKDE_soil_runoff()
        py"""
        import numpy as np

        def LoadNumpyArray():
            KDE_runoff = np.load(os.getcwd() + "/KDE/KDE_runoff_soil_complete_3D.npy")
            runoff_vector = np.load(os.getcwd() + "/KDE/runoff_range_3D_complete.npy")
            soil_vector = np.load(os.getcwd() + "/KDE/soil_range_3D_complete_runoff.npy")
            z_integrated = np.sum(KDE_runoff[:,:], axis = 0)
            z_norm = KDE_runoff[:,:]/z_integrated[None, :]
            print(np.sum(z_norm[:,:], axis = 0))
            return z_norm, soil_vector, runoff_vector
        """
        z_norm, soil_vector, runoff_vector = py"LoadNumpyArray"()
    end
end
