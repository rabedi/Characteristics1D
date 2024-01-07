classdef PF
    properties

        % inputs

        % normalized loading rate: if < 0, values are plotted versus y the
        % master-curve for any loading rate
        ap = -1;
        %CZM_normalization4AT1_2: means that when we use AT1/AT2 models, we change
        %the normalization from sigmaNormalization = sqrt(G/Eb) and
        %lengthNormalization = b to sigmaNormalization = sigma_coh, and
        %lengthNormalization = l_coh. This option = 1 is useful if we compare AT1
        %and AT2 with CZM models
        CZM_normalization4AT1_2 = 0;
        
        % for P use 0, for H use 1, for E use -1
        isHyper = 1;
        % AT2: xi = 0,  AT1: xi = 1, CZM: xi = 2
        xi = 1;
        % omega = (1 - D)^2 for omegaCZM == 0, else is the rational
        % expression
        % if -1, CZM model is decided based on xi
        omegaCZM = -1;

        % for omegaCZM == 1
        CZM_model_name = 'Linear'; %'Bilinear', 'Exponential', 'Hyperbolic', 'Concrete'

        bPrime = -1; % the right value is chosen for each model if bPrime < 0, otherwise, the given value is used
        % the analytical solution is only valid for parabolic solution. The
        % gamma-based solution for the hyperbolic case is incorrect
        useAnalyticalSln = 0;


        % computed values
        % vectors
        epsilon_p;
        sigma_p;
        D_vec;
        Dp_vec;
        Dpp_vec;
        omegaD_vec;

        sigmap_Max, sigmap_Max_eps, eps_f, phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain;

        y_vecs; % contains sigma_p, ....
        scalars; % contains sigmaP_Max, ....
        vecs_names = {'sigmap', 'D', 'Dp', 'Dpp', 'omega'};
        scalar_names = {'sigma_m', 'epsilon_y_m', 'epsilon_y_f', 'phi', 'phi_unloading', 'phi_loading', 'brittleness_phi', 'brittleness_strain'};
    end
    methods
        function objout = Compute(obj)
            objout = obj;
            if (objout.omegaCZM < 0)
                if (objout.xi == 2)
                    objout.omegaCZM = 1;
                else
                    objout.omegaCZM = 0;
                end
            end
            if (objout.isHyper >= 0) % dynamic models
                if (objout.ap > 0)
                    [objout.epsilon_p, objout.sigma_p, objout.D_vec, objout.Dp_vec, objout.Dpp_vec, objout.omegaD_vec, objout.sigmap_Max, objout.sigmap_Max_eps, objout.eps_f, ...
                    objout.phi, objout.phi_unloading, objout.phi_loading, objout.brittleness_phi, objout.brittleness_strain] = ...
                    ComputePF_withLoadingRate(objout.ap, objout.useAnalyticalSln, objout.isHyper, objout.xi, objout.omegaCZM, objout.bPrime, objout.CZM_model_name, objout.CZM_normalization4AT1_2);
                else
                    [objout.epsilon_p, objout.sigma_p, objout.D_vec, objout.Dp_vec, objout.Dpp_vec, objout.omegaD_vec, objout.sigmap_Max, objout.sigmap_Max_eps, objout.eps_f, ...
                    objout.phi, objout.phi_unloading, objout.phi_loading, objout.brittleness_phi, objout.brittleness_strain] = ...
                    ComputePF(objout.useAnalyticalSln, objout.isHyper, objout.xi, objout.omegaCZM, objout.bPrime, objout.CZM_model_name, objout.CZM_normalization4AT1_2);
                end
            else
                    [objout.epsilon_p, objout.sigma_p, objout.D_vec, objout.Dp_vec, objout.Dpp_vec, objout.omegaD_vec, objout.sigmap_Max, objout.sigmap_Max_eps, objout.eps_f, ...
                    objout.phi, objout.phi_unloading, objout.phi_loading, objout.brittleness_phi, objout.brittleness_strain] = ...
                    ComputePF_Elliptic(objout.xi, objout.omegaCZM, objout.bPrime, objout.CZM_model_name, objout.CZM_normalization4AT1_2);
            end
            objout.y_vecs{1} = objout.sigma_p;
            objout.y_vecs{2} = objout.D_vec;
            objout.y_vecs{3} = objout.Dp_vec;
            objout.y_vecs{4} = objout.Dpp_vec;
            objout.y_vecs{5} = objout.omegaD_vec;
            objout.scalars = [objout.sigmap_Max, objout.sigmap_Max_eps, objout.eps_f, objout.phi, objout.phi_unloading, objout.phi_loading, objout.brittleness_phi, objout.brittleness_strain];
        end
    end
end