classdef PF
    properties
        % this is introduced to resolve the problem with shifting eps_final
        % behavior in convergence runs, -1 means it's used to prevent a shift in Dfinal 
        DMaxHyper = -1;

        % inputs
        % asymptoticMode = 
        %   2  ->    high loading rate approximation AND only the asymptotic part
        %   1  ->    high loading rate approximation
        %   -1 ->   elliptic limit
        %   0  ->   compute the exact solution
        asymptoticMode = 1;

        % damping factor w.r.t. time scale used
        dampingFactor = 1.0;
        % k parameter
        waveSpeedD_to_waveSpeedRef = 1.0;
        % log of wave speed ratio, if > -100, wave speed ratio is computed
        % from this
        l_waveSpeedD_to_waveSpeedRef = -101;

        % normalized loading rate: if < 0, values are plotted versus y the
        % master-curve for any loading rate
        lap = -101; % if lap < -100, it's not used, if > -100, ap is computed from that
        ap = -1;
        %CZM_normalization4AT1_2: means that when we use AT1/AT2 models, we change
        %the normalization from sigmaNormalization = sqrt(G/Eb) and
        %lengthNormalization = b to sigmaNormalization = sigma_coh, and
        %lengthNormalization = l_coh. This option = 1 is useful if we compare AT1
        %and AT2 with CZM models
        CZM_normalization4AT1_2 = 1;
        
        % for P use 0, for H use 1, for E use -1
        isHyper = 1;
        % AT2: xi = 0,  AT1: xi = 1, CZM: xi = 2
        xi = 1;
        % omega = (1 - D)^2 for omegaCZM == 0, else is the rational
        % expression
        % if -1, CZM model is decided based on xi
        omegaCZM = -1;

% tauModel = 
%           auto            -> if CZM_normalization4AT = 1 ->   tauModel = 'lcoh2cd'
%                              else                             tauModel = 'b2cd'

%           b2cd            tauScale = b/cd, cd = Phase field speed
%           lcoh2cd         tauScale = lcoh2cd/cd, cd = Phase field speed
        tauModel = 'b2cd'; %'auto';

        % for omegaCZM == 1
        CZM_model_name = 'Linear'; %'Bilinear', 'Exponential', 'Hyperbolic', 'Concrete'

        bPrime = -1; % the right value is chosen for each model if bPrime < 0, otherwise, the given value is used
        % the analytical solution is only valid for parabolic solution. The
        % gamma-based solution for the hyperbolic case is incorrect
        bTimesPiCZM = 1.0; % for CZM-Wu value of 2 is max (sharp drop), 1 is good, 0.5 is two times slower decay, etc.

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
        vecs_names = {'sigmap', 'D', 'Dp', 'Dpp', 'omega', 'epsilonp'};
        vecs_names_latex = {'$$ \bar{\sigma} $$', '$$ d $$', '$$ d^\prime $$', '$$ d^{\prime\prime} $$', '$$ \omega(d) $$', '$$\bar{\epsilon} $$'};
        scalar_names = {'sigma_m', 'epsilon_y_m', 'epsilon_y_f', 'phi', 'phi_unloading', 'phi_loading', 'brittleness_phi', 'brittleness_strain'};
        scalar_names_latex = {'\bar{\sigma}_M', '\bar{\epsilon}_M', '\bar{\epsilon}_f', ...
            '\bar{\phi}', '\bar{\phi}_u', '\bar{\phi}_M', 'B^\phi', 'B^\epsilon'};
        % a's for CZM-based omega
        pis2 = 1;
        p = 2.0;
        a1 = 0.0, a2 = 0.0, a3 = 0.0;

        cAlpha = 0;

        % 1 - D < tol -> considered full damage
        delD = 1e-3;
        % time step of epsilon
        del_eps = 1e-4;
        hasElasticRegime;
        E = 1.0; % normalized elastic modulus
    end
    methods
        function objout = Compute(obj, plotResults)
            if nargin < 2
                plotResults = 0;
            end
            objout = obj.ComputeNew(plotResults);
%            objout = obj.Compute_Old_version1();
        end
        function objout = ComputeNew(obj, plotResults)
            if nargin < 2
                plotResults = 0;
            end
            % ellptic cases
            if (obj.isHyper < 0)
                obj.asymptoticMode = -1;
            end
            if (obj.asymptoticMode == -1)
                obj.isHyper = -1;
            end

            % dynamic solution with ap < 0 -> only asmyptotic solution is available
            if (obj.isHyper >= 0) 
                if (obj.ap < 0)
                    obj.asymptoticMode = 2;
                end
            end

            asymptoticModeTmp = obj.asymptoticMode;
            if (obj.omegaCZM < 0)
                if (obj.xi == 2)
                    obj.omegaCZM = 1;
                else
                    obj.omegaCZM = 0;
                end
            end

            epsDot = obj.ap;
            onlyKeepAsymptotic = (asymptoticModeTmp == 2);
            if (onlyKeepAsymptotic)
                asymptoticModeTmp = 1;
                epsDot = 1.0;
            end
            objout = Initialize_Stage2(obj);
            wspdInv = 1.0 / objout.waveSpeedD_to_waveSpeedRef;
            bbarCorrection = 0;
            % if tauModel == lcoh2cd, it means that the length scale =
            % lcoh, so tau_d / tauScale = (b/c_d) / (lcoh / c_u) = 
            % (b/lcoh) / (c_d / c_u) = bPrime / k
            % in this case wspdInv must multiply wspdInv

            if (~onlyKeepAsymptotic)
                bbarCorrection = (strcmp(objout.tauModel, 'lcoh2cd') == 1);
                if (strcmp(objout.tauModel, 'auto') == 1)
                    bbarCorrection = objout.CZM_normalization4AT1_2;
                end
            end
            if (bbarCorrection)
                wspdInv = wspdInv * objout.bPrime;
            end

            dDotCoef = objout.dampingFactor * wspdInv;
            dDDotCoef = wspdInv * wspdInv;

            if ((obj.isHyper == -1) || (asymptoticModeTmp == -1)) % elliptic OR aymptoticMode elliptic
                       [epsilon_p_vec, sigma_p_vec, D_vec, Dp_vec, Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_eps, eps_f, ...
                phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
                ComputePF_Elliptic(objout.xi, objout.omegaCZM, objout.bPrime, objout.CZM_model_name, objout.CZM_normalization4AT1_2);
                objout.epsilon_p = epsilon_p_vec;
                objout.sigma_p = sigma_p_vec;
                objout.D_vec = D_vec;
                objout.Dp_vec = Dp_vec;
                objout.Dpp_vec = Dpp_vec;
                objout.omegaD_vec = omegaD_vec;
            else
                eps0 = 0.0;
                if (obj.isHyper == 0) % parabolic
                    [objout, eps0] = objout.SolveParabolic(epsDot, dDotCoef, asymptoticModeTmp);
                elseif (obj.isHyper == 1) % hyperbolic
                    [objout, eps0] = objout.SolveHyperbolic(epsDot, dDotCoef, dDDotCoef, asymptoticModeTmp);
                end
                if ((onlyKeepAsymptotic) && (objout.hasElasticRegime)) % remove the first point if needed and update the stres
                    delSigma = -eps0 * objout.omegaD_vec;
                    objout.sigma_p = objout.sigma_p + delSigma;
                    objout.epsilon_p = objout.epsilon_p  - eps0;
                    sz = length(objout.epsilon_p);
                    objout.epsilon_p = objout.epsilon_p(2:sz);
                    objout.sigma_p = objout.sigma_p(2:sz);
                    objout.D_vec = objout.D_vec(2:sz);
                    objout.Dp_vec = objout.Dp_vec(2:sz);
                    objout.Dpp_vec = objout.Dpp_vec(2:sz);
                    objout.omegaD_vec = objout.omegaD_vec(2:sz);
                end
            end

            % compute max stress and get energies
            sz = length(objout.sigma_p);
            [objout.sigmap_Max, i] = max(objout.sigma_p);
            objout.sigmap_Max_eps = objout.epsilon_p(i);
            objout.eps_f = objout.epsilon_p(sz);
            objout.phi = trapz(objout.epsilon_p, objout.sigma_p);
            objout.phi_unloading = trapz(objout.epsilon_p(i:sz), objout.sigma_p(i:sz));
            objout.phi_loading = objout.phi - objout.phi_unloading;
            objout.brittleness_phi = objout.phi_loading / objout.phi;
            objout.brittleness_strain = objout.sigmap_Max_eps / objout.eps_f;

            objout.y_vecs{1} = objout.sigma_p;
            objout.y_vecs{2} = objout.D_vec;
            objout.y_vecs{3} = objout.Dp_vec;
            objout.y_vecs{4} = objout.Dpp_vec;
            objout.y_vecs{5} = objout.omegaD_vec;
            objout.y_vecs{6} = objout.epsilon_p;
            objout.scalars = [objout.sigmap_Max, objout.sigmap_Max_eps, objout.eps_f, objout.phi, objout.phi_unloading, objout.phi_loading, objout.brittleness_phi, objout.brittleness_strain];

            if (plotResults)
                sc = objout.scalars
                epss = objout.epsilon_p;
                sigs = objout.sigma_p;
                Ds = objout.D_vec;
    
                figure(1);
                plot(epss, sigs);
                figure(2);
                plot(epss, Ds);
            end
        end
        function [objout, str_vec_out] = Initialize_Stage1(obj, vec_In)
            [objout, str_vec_out] = obj.Initialize_Stage1_Aux(vec_In{1}, vec_In{2}, vec_In{3}, vec_In{4}, vec_In{5}, vec_In{6}, vec_In{7}, vec_In{8}, vec_In{9}, vec_In{10});
        end
        function [objout, str_vec_out] = Initialize_Stage1_Aux(obj, isHyperIn, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4ATIn, bTimesPiCZM_s, CZM_modelNameIn, tauModelIn)
            if (strcmp(la_s, 'none') == 0)
                lav = str2num(la_s);
                if ((length(lav) == 0) || (lav < -99))
                    isHyperIn = -1;
                end
            end
            obj.tauModel = tauModelIn;
            obj.isHyper = isHyperIn;
            obj.lap = -101;
            obj.ap = -1;
            la = nan;
            l_cDRatio = nan;
            df = nan;
            bTimesPiCZM = nan;
            CZM_modelName = 'none';
            CZM_normalization4AT = CZM_normalization4ATIn;
            if (isHyperIn == -1)
                obj.asymptoticMode = -1;
                isApproximate = 0;
                l_cD2c = nan;
                df = nan;
            else
                la = str2num(la_s);
                hasLoadingRate = (length(la) > 0);
                if (la <= -100)
                    hasLoadingRate = 0;
                end
                if (hasLoadingRate)
                    obj.lap = la;
                    obj.ap = power(10, la);
                    if (isApproximate == 1)
                        obj.asymptoticMode = 1;
                    else
                        obj.asymptoticMode = 0;
                    end

                    l_cD2c = str2num(l_cD2c_s);
                    hascDRatio = (length(l_cD2c) > 0);
                    if (l_cD2c <= -100)
                        hascDRatio = 0;
                    end
                    if (hascDRatio)
                        l_cDRatio = l_cD2c;
                        obj.l_waveSpeedD_to_waveSpeedRef = l_cD2c;
                    end
                    df_v = str2num(df_s);
                    if (df_v > 0)
                        df = df_v;
                        obj.dampingFactor = df_v;
                    end
                    
                else
                    isApproximate = 1;
                    obj.asymptoticMode = 2;
                end
            end
            if (contains(model_s, 'CZM'))
                CZM_normalization4AT = 1;
                obj.omegaCZM = 1;
                if (contains(model_s, 'W'))
                    obj.xi = 2;
                elseif (contains(model_s, 'L'))
                    obj.xi = 1;
                end
                bTimesPiCZM_v = str2num(bTimesPiCZM_s);
                if (length(bTimesPiCZM_v) > 0)
                    bTimesPiCZM = bTimesPiCZM_v;
                    obj.bTimesPiCZM = bTimesPiCZM_v;
                end
                CZM_modelName = 'Linear';
                if (strcmp(CZM_modelNameIn, 'Bilinear') == 1)
                    CZM_modelName = 'Bilinear';
                elseif (strcmp(CZM_modelNameIn, 'Exponential') == 1)
                    CZM_modelName = 'Exponential';
                elseif (strcmp(CZM_modelNameIn, 'Hyperbolic') == 1)
                    CZM_modelName = 'Hyperbolic';
                elseif (strcmp(CZM_modelNameIn, 'Concrete') == 1)
                    CZM_modelName = 'Concrete';
                end
                obj.CZM_model_name = CZM_modelName;
            else
                obj.omegaCZM = 0;
                if (strcmp(model_s, 'AT1') == 1)
                    obj.xi = 1;
                elseif (strcmp(model_s, 'AT2') == 1)
                    obj.xi = 0;
                end
                if (CZM_normalization4AT == 0)
                    obj.CZM_normalization4AT1_2 = 0;
                else
                    obj.CZM_normalization4AT1_2 = 1;
                end
            end
            cntr = 1;
            str_vec_out{cntr} = num2str(isHyperIn);
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(isApproximate);
            cntr = cntr + 1;
            str_vec_out{cntr} = model_s;
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(la);
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(l_cDRatio);
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(df_s);
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(CZM_normalization4AT);
            cntr = cntr + 1;
            str_vec_out{cntr} = num2str(bTimesPiCZM);
            cntr = cntr + 1;
            str_vec_out{cntr} = CZM_modelName;
            cntr = cntr + 1;
            str_vec_out{cntr} = obj.tauModel;
            objout = obj;
        end

        function objout = Initialize_Stage2(obj)
            obj.hasElasticRegime = 1;
            if (obj.lap > -100)
                obj.ap = power(10, obj.lap);
            end
            if (obj.l_waveSpeedD_to_waveSpeedRef > -100)
                obj.waveSpeedD_to_waveSpeedRef = power(10.0, obj.l_waveSpeedD_to_waveSpeedRef);
            end
            if (obj.omegaCZM <= 0) 
                if (obj.xi == 0) % AT2
                    obj.hasElasticRegime = 0;
                    obj.cAlpha = 2.0;
                else
                    obj.cAlpha = 8.0 / 3.0;
                end

                if (obj.bPrime < 0)
                    if (obj.CZM_normalization4AT1_2 == 1)
                        if (obj.xi == 0) % AT2
                            obj.bPrime = 27.0 / 256.0;
                        elseif (obj.xi == 1) % AT1
                            obj.bPrime = 3.0 / 8.0;
                        end
                    else
                        if (obj.bPrime < 0)
                            obj.bPrime = 1.0;
                        end
                    end
                end
            else % CZM's
                if (obj.xi == 2)
                    if (obj.bTimesPiCZM > 0)
                        obj.bPrime = 1.0 / pi * obj.bTimesPiCZM;
                    end
                    obj.cAlpha = pi;
                    obj.a1 = 4.0 / pi / obj.bPrime; 
                elseif (obj.xi == 1) % Lorentz
                    if (obj.bPrime < 0)
                        if (obj.bTimesPiCZM > 0)
                            obj.bPrime = 0.25 * obj.bTimesPiCZM;
                        else
                            obj.bPrime = 1.0;
                        end
                    end
                    obj.cAlpha = 8.0 / 3.0;
                    obj.a1 = 3.0 / 4.0 / obj.bPrime; 
                end
            end

            CZM_model_name = obj.CZM_model_name;
            p = 2;
            pis2 = 1;
            a2 = 0.0;
            a3 = 0.0;
            if (obj.xi == 2)
                if (strcmp(CZM_model_name, 'Linear') == 1)
                    a2 = -0.5;
                elseif (strcmp(CZM_model_name, 'Bilinear') == 1)
                    a2 = 0.03687;
                    a3 = 20.8343;
                elseif (strcmp(CZM_model_name, 'Exponential') == 1)
                    p = 2.5;    
                    pis2 = 0;
                    a2 = power(2.0, 5.0/3.0) - 3.0;
                elseif (strcmp(CZM_model_name, 'Hyperbolic') == 1)
                    p = 4;
                    pis2 = 0;
                    a2 = power(2.0, 7.0/3.0) - 4.5;
                elseif (strcmp(CZM_model_name, 'Concrete') == 1)
                    a2 = 1.3868;
                    a3 = 0.6567;
                end
            elseif (obj.xi == 1)
                    a2 = 1.0;
            end

            obj.pis2 = pis2;
            obj.p = p;
            obj.a2 = a2;
            obj.a3 = a3;

            objout = obj;
        end



        function [omegaD, omegaPD, half_alphaP] = ComputeOmegas_AlphaP(obj, D)
            if (1 - D <  0)
                omegaD = 0.0;
                omegaPD = 1.0;
                half_alphaP = 1.0;
            end

            if (obj.omegaCZM <= 0)
                omegaD = (1.0 - D);
                omegaPD = -2.0 * omegaD;
                omegaD = omegaD  * omegaD;
                if (obj.xi == 0) % AT2
                    half_alphaP = D;
                elseif (obj.xi == 1) % AT1
                    half_alphaP = 0.5;
                end
                return;
            end
            factor = obj.a1;
%            D = 0.6;
            if (obj.xi == 2) % CZM-Wu
                half_alphaP = 1 - D;
            else
                half_alphaP = 0.5;
            end
            omd_power = 1.0 - D;
            if (obj.pis2 == 1)
                omd_power_p = -2.0 * omd_power;
                omd_power = omd_power * omd_power;
            else
                omd_power_p = power(omd_power, obj.p - 1);
                omd_power = omd_power * omd_power_p; 
                omd_power_p = -obj.p * omd_power_p;
            end
            Pd = 1.0 + obj.a2 * D * (1.0 + obj.a3 * D);
            Pdp = obj.a2 * (1.0 + 2.0 * D * obj.a3);
            denom = (omd_power + factor * D * Pd);
            omegaD = omd_power / denom;
            num = factor * (D * (omd_power_p * Pd - omd_power * Pdp) - omd_power * Pd);
            omegaPD = num / denom / denom;
        end
        % use_HLR: high loading rate appeoximation
        function [rhs, omegaD, omegaPD] = Compute_PF_RHS_omega_alphaPrime(obj, D, eps, eps_ini, use_HLR)
            [omegaD, omegaPD, half_alphaP] = ComputeOmegas_AlphaP(obj, D);
            if (~use_HLR)
                rhs = -0.25 * obj.cAlpha * omegaPD * obj.bPrime * obj.E * eps * eps - half_alphaP;
            else
                eps = eps - eps_ini;
                rhs = -0.25 * obj.cAlpha * omegaPD * obj.bPrime * obj.E * eps * eps;
            end
            if (rhs < 0.0)
                rhs = 0.0;
            end
        end

        function [objout, eps0] = SolveParabolic(obj, epsDot, dDotCoef, asymptoticModeTmp)
            objout = obj;
            Dmax = 1.0 - obj.delD;
%            Dmax = 1 - 0.01;
            deleps = obj.del_eps;

            use_HLR = (asymptoticModeTmp == 1);

            if (deleps < 0)
                if (epsDot > 10)
                    deleps = -deleps * power(epsDot, 2.0/3.0);
                else
                    deleps = -deleps;
                end
            else
                if (epsDot > 10000)
                    deleps = max(2e-3, deleps);
                elseif (epsDot > 1000)
                    deleps = max(3e-4, deleps);
                end
            end
            
            % equation is cast in terms of epsilon
            dEpsCoefInv = 1.0 / (dDotCoef * epsDot);

            epss = 0;
            sigs = 0;
            Ds = 0;
            Dps = 0;
            Dpps = 0;
            omegaDs = 1;
            cntr = 1;
            
            % calculating initiation value
            eps0 = 0.0;
            if (obj.hasElasticRegime)
                [omegaD0, omegaPD0, half_alphaP0] = obj.ComputeOmegas_AlphaP(0.0);
                eps0 = sqrt(-4.0 * half_alphaP0 / obj.cAlpha / obj.bPrime / omegaPD0);
                cntr = 2;
                epss(cntr) = eps0;
                sigs(cntr) = obj.E * eps0;
                Ds(cntr) = 0;
                Dps(cntr) = 0;
                Dpps(cntr) = 0;
                omegaDs(cntr) = omegaD0;
            end

            delx = deleps;
            half_delx = 0.5 * delx;
            delx_div6 = delx / 6.0;
            delxInv = 1.0 / delx;
            delxInv2 = delxInv * delxInv;

            Dg1 = 0;
            while ((Ds(cntr) < Dmax) && (cntr < 160000000))
                x = epss(cntr);
                y = Ds(cntr);

                pr = 0;
                if (0)
%                if ((y >= 0.1) && (Dg1 == 0))
                    fprintf(1, 'start debugging\n');
                    pr = 1;
                    Dg1 = 1;
                end
                xNew = x + delx;

                % applying RK4:
    
                % y = D, 
                % x = epsilon;
                % dy/dx = f(x, y) 
                % y(x + delx) = y(x) + delx/6 * (k1 + 2k2 + 2k3 + k4)
                % k1 = f(x, y)
                % k2 = f(x + delx / 2, y + delx / 2 * k1)
                % k3 = f(x + delx / 2, y + delx / 2 * k2)
                % k4 = f(x + delx    , y + delx k3    )

                % k1
                [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y, x, eps0, use_HLR);
                k1 = rhs * dEpsCoefInv; 
                if (pr == 1)
                    D = y
                    eps = x
                    bPrimePi = obj.bPrime * pi
                    omegaD
                    omegaPD
                    tmpD = omegaPD * obj.bPrime
                    rhs
                    sigma2 = obj.E * omegaD * eps
                    a = 12;
                end


                % calculating k2
                x2Use = x + half_delx;
                y2Use = y + half_delx * k1;
                DisLessThan1 = (y2Use <= Dmax);
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use, x2Use, eps0, use_HLR);
                    k2 = rhs * dEpsCoefInv; 
                    % x2Use = x + half_delx;
                    y2Use = y + half_delx * k2;
                    DisLessThan1 = (y2Use <= Dmax);
                end

                % calculating k3
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use, x2Use, eps0, use_HLR);
                    k3 = rhs * dEpsCoefInv; 
                    x2Use = x + delx;
                    y2Use = y + delx * k3;
                    DisLessThan1 = (y2Use <= Dmax);
                end

                % calculating k4
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use, x2Use, eps0, use_HLR);
                    k4 = rhs * dEpsCoefInv; 
                    yNew = y + delx_div6 * (k1 + k4 + 2.0 * (k2 + k3));
                    DisLessThan1 = (yNew <= Dmax);
                else
                    yNew = 1.0;
                    sigma = 0.0;
                end
                eps = xNew;
                D = 1.0;
                sigma = 0.0;
                omegaD = 0.0;
                if (DisLessThan1)
                    D = yNew;
                    [omegaD, omegaPD, half_alphaP] = ComputeOmegas_AlphaP(obj, D);
                    sigma = obj.E * omegaD * eps;
                end
                
                cntr = cntr + 1;

                epss(cntr) = eps;
                sigs(cntr) = sigma;
                Ds(cntr) = D;
                Dprev = Ds(cntr - 1);
                Dps(cntr) = delxInv * (D - Dprev);
                Dpp = 0.0;
                if (cntr > 3)
                    % this is technically Dpp of last step ...
                    Dpp = (D + Ds(cntr - 2) - 2.0 * Dprev) * delxInv2;
                end
                Dpps(cntr) = Dpp;
                omegaDs(cntr) = omegaD;
            end
            sz = length(epss);
            Dps(sz) = Dps(sz - 1);
            Dpps(sz - 1) = Dpps(sz - 2);
            Dpps(sz) = Dpps(sz - 2);
            objout.epsilon_p = epss;
            objout.sigma_p = sigs;
            objout.D_vec = Ds;
            objout.Dp_vec = Dps;
            objout.Dpp_vec = Dpps;
            objout.omegaD_vec = omegaDs;
        end


        function [objout, eps0] = SolveHyperbolic(obj, epsDot, dDotCoef, dDDotCoef, asymptoticModeTmp)
            iterMax = 8000000;
            deleps = obj.del_eps;
            if (deleps < 0)
                if (epsDot > 10)
                    deleps = -deleps * epsDot;
                else
                    deleps = -deleps;
                end
            else
                if (epsDot > 10000)
                    deleps = max(2e-3, deleps);
                elseif (epsDot > 1000)
                    deleps = max(3e-4, deleps);
                end
            end
            % baseline deleps
            success = 0;
            maxTry = 10;
            tryNo = 0;
            while ((success == 0) && (tryNo < maxTry))
                [objout, eps0, badD, Dmax] = obj.SolveHyperbolicAux(epsDot, dDotCoef, dDDotCoef, asymptoticModeTmp, deleps, iterMax);
                if (badD)
                    deleps = deleps * 0.1;
                    iterMax = iterMax * 10;
                else
                    Dmx = max(objout.D_vec);
                    DmxLim = 0.99;
                    if (obj.DMaxHyper > 0)
                        DmxLim = obj.DMaxHyper;
                    end
                    DmxLim = 0.999999 * DmxLim;
                    success = (Dmx > DmxLim);
                    if (tryNo > 5)
                            deleps = deleps * 5;
                    end
                    if (success == 0)
                        iterMax = iterMax * 10;
                        tp1 = tryNo + 1;
                        fprintf(1, 'bPrime = %g, xi = %d, tryNo = %d, iterMax = %d, Dmax %g\n', obj.bPrime, obj.xi, tp1, iterMax, Dmax);
                    end
                end
                tryNo = tryNo + 1;
            end
            if (success == 0)
                fprintf(1, 'Cannot solve this set-up\n');
                pause;
            end
        end
        function [objout, eps0, badD, Dmax] = SolveHyperbolicAux(obj, epsDot, dDotCoef, dDDotCoef, asymptoticModeTmp, deleps, iterMax)

            badD = 0;
            objout = obj;
            Dmax = 1.0;
            % not needed to stop before 1 as Dmax reaches 1
%            Dmax = 1.0 - obj.delD;

            use_HLR = (asymptoticModeTmp == 1);


            % equation is cast in terms of epsilon
            % C2 d2 D / d eps2 + C1 d D / d esp
            C2 = dDDotCoef * epsDot * epsDot;
            C1 = dDotCoef * epsDot;
            if (use_HLR) % for asymptotic solution, we no longer have the damping term
                C1 = 0.0;
            end
            C2_inv = 1.0 / C2;
            C1_div_C2 = C1 * C2_inv;

            sC2 = sqrt(C2);
            if (obj.asymptoticMode ~= 2) ...
                    if ((sC2 <= 0.11) || (C1 > 6 * C2))
                        Dmax = 1.0 - obj.delD;
                    elseif ((sC2 <= 0.41) || (C1 > 4 * C2))
                        Dmax = 1.0 - 0.3 * obj.delD;
                    end
            end
            if (obj.DMaxHyper > 0)
                Dmax = obj.DMaxHyper;
            end
            % C2 d2 D / d eps2 + C1 d D / d eps = rhs ->

            % C2 d Dp / d eps = rhs - C1 Dp
            % d D / d eps = Dp
            
            % -> 

            % d Dp / d eps = 1/C2 (rhs - C1 Dp) = -C1/C2 Dp + 1/C2 rhs
            % d D / d eps = Dp

            epss = 0;
            sigs = 0;
            Ds = 0;
            Dps = 0;
            Dpps = 0;
            omegaDs = 1;
            cntr = 1;
            
            % calculating initiation value
            eps0 = 0.0;
            if (obj.hasElasticRegime)
                [omegaD0, omegaPD0, half_alphaP0] = obj.ComputeOmegas_AlphaP(0.0);
                eps0 = sqrt(-4.0 * half_alphaP0 / obj.cAlpha / obj.bPrime / omegaPD0);
                cntr = 2;
                epss(cntr) = eps0;
                sigs(cntr) = obj.E * eps0;
                Ds(cntr) = 0;
                Dps(cntr) = 0;
                Dpps(cntr) = 0;
                omegaDs(cntr) = omegaD0;
            end

            delx = deleps;
            half_delx = 0.5 * delx;
            delx_div6 = delx / 6.0;
            delxInv = 1.0 / delx;
            delxInv2 = delxInv * delxInv;

            while ((Ds(cntr) < Dmax) && (cntr < iterMax))
                x = epss(cntr);
                y(1) = Dps(cntr);
                y(2) = Ds(cntr);
                y2Use = y;

                xNew = x + delx;

                % applying RK4:
    
                % y = D, 
                % x = epsilon;
                % dy/dx = f(x, y) 
                % y(x + delx) = y(x) + delx/6 * (k1 + 2k2 + 2k3 + k4)
                % k1 = f(x, y)
                % k2 = f(x + delx / 2, y + delx / 2 * k1)
                % k3 = f(x + delx / 2, y + delx / 2 * k2)
                % k4 = f(x + delx    , y + delx k3    )

                % k1
                [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use(2), x, eps0, use_HLR);
                k1(1) = -C1_div_C2 * y2Use(1) + C2_inv * rhs;
                k1(2) = y2Use(1);

                % calculating k2
                x2Use = x + half_delx;
                y2Use = y + half_delx * k1;
                DisLessThan1 = (y2Use(2) <= Dmax);
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use(2), x2Use, eps0, use_HLR);
                    k2(1) = -C1_div_C2 * y2Use(1) + C2_inv * rhs;
                    k2(2) = y2Use(1);
                    % x2Use = x + half_delx;
                    y2Use = y + half_delx * k2;
                    DisLessThan1 = (y2Use(2) <= Dmax);
                end

                % calculating k3
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use(2), x2Use, eps0, use_HLR);
                    k3(1) = -C1_div_C2 * y2Use(1) + C2_inv * rhs;
                    k3(2) = y2Use(1);
                    x2Use = x + delx;
                    y2Use = y + delx * k3;
                    DisLessThan1 = (y2Use(2) <= Dmax);
                end

                % calculating k4
                if (DisLessThan1)
                    [rhs, omegaD, omegaPD] = obj.Compute_PF_RHS_omega_alphaPrime(y2Use(2), x2Use, eps0, use_HLR);
                    k4(1) = -C1_div_C2 * y2Use(1) + C2_inv * rhs;
                    k4(2) = y2Use(1);
                    yNew = y + delx_div6 * (k1 + k4 + 2.0 * (k2 + k3));
                    Dtmp = yNew(2);
                    DisLessThan1 = (yNew(2) <= Dmax);
                else
                    yNew(1) = 0.0;
                    yNew(2) = 1.0;
                    sigma = 0.0;
                end
                eps = xNew;
                D = 1.0;
                Dp = 0.0;
                sigma = 0.0;
                omegaD = 0.0;
                if (DisLessThan1)
                    D = yNew(2);
                    Dp = yNew(1);
                    [omegaD, omegaPD, half_alphaP] = ComputeOmegas_AlphaP(obj, D);
                    sigma = obj.E * omegaD * eps;
                end
                
                cntr = cntr + 1;

                epss(cntr) = eps;
                sigs(cntr) = sigma;
                if ((D < 0.0) || ((cntr > 1) && (D < 0.999 * Ds(cntr - 1))))
                    badD = 1;
                    return;
                end
                Ds(cntr) = D;
                Dps(cntr) = Dp;
                Dpprev = Dps(cntr - 1);
                Dpps(cntr) = delxInv * (Dp - Dpprev);
                omegaDs(cntr) = omegaD;
            end

            sz = length(epss);
            Dps(sz) = Dps(sz - 1);
            Dpps(sz - 1) = Dpps(sz - 2);
            Dpps(sz) = Dpps(sz - 2);
            
            objout.epsilon_p = epss;
            objout.sigma_p = sigs;
            objout.D_vec = Ds;
            objout.Dp_vec = Dps;
            objout.Dpp_vec = Dpps;
            objout.omegaD_vec = omegaDs;
        end
        % does not support real solution [asymptoticMode == 0 is not
        % supported]
        function objout = Compute_Old_version1(obj)
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