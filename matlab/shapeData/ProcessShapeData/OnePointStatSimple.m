classdef OnePointStatSimple
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% raw data
        %%%% y_i(x_i) % y is vals, x is locs
        vals;  
        num = -1; % actual size of vals
        % point values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% basic statistics
        minV = nan;
        minLoc = nan;

        maxV = nan;
        maxLoc = nan;
                
        % use angleOption to compute the following
        meanV = nan;
        sdivV = nan;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% derived statistics
        covV; % nsdiv
        spanV; % max - min
        nspanV; 
        pdf_xAxis = [];
        pdf_yAxis = [];
        cdf_xAxis = [];
        cdf_yAxis = [];

        %%%%%%%%%%%%%%%%
        Weibull_shape = 0;
        Weibull_scale = 0;
        Weibull_N = 0;
        Weibull_h = 0;
        Weibull_p = 0;
        
        % my calculated vals
        WeibullR_shape = 0;
        WeibullR_scale = 0;
        WeibullR_N = 0;
        WeibullR_R2 = 0;
        WeibullR_X = 0;
        WeibullR_Y = 0;
        WeibullR_xreg = 0;
        WeibullR_yreg = 0;
        WeibullR_p = 0; % [slope, intercept]

        Gauss_mu = 0;
        Gauss_sigma = 0;
        Gauss_h = 0;
        Gauss_p = 0;
    end
    methods
        function objout = ComputeVals(obj, valsIn, min4Weibull, doWeibull, doGauss)
            if (nargin < 2)
                valsIn = [];
            end
            if (nargin < 3)
                min4Weibull = 0.5;
            end
            if (nargin < 4)
                doWeibull = 0;
            end
            if (nargin < 5)
                doGauss = 0;
            end
            if (isempty(obj.vals))
                obj.vals = valsIn;
            end
            obj.num = length(obj.vals);
            [obj.minV, obj.minLoc] = min(obj.vals);
            [obj.maxV, obj.maxLoc] = max(obj.vals);
            obj.meanV = mean(obj.vals);
            obj.sdivV = std(obj.vals);
            obj.covV = obj.sdivV * obj.sdivV;
            obj.spanV = obj.maxV - obj.minV;
            obj.nspanV = obj.spanV / obj.meanV;
            [obj.pdf_yAxis, obj.pdf_xAxis] = ksdensity(obj.vals);
            [obj.cdf_yAxis, obj.cdf_xAxis] = ksdensity(obj.vals, 'Function','cdf');
            % Weibull fit
            if (doWeibull)
                data = obj.vals - min4Weibull + 1e-9;
                minV = min(data);
                if (minV <= 0)
                    min4Weibull
                    minV
                    pause;
                end
                [paramEsts, paramCIs] = wblfit(data);   
                obj.Weibull_shape = paramEsts(2); % (1)
                obj.Weibull_scale = paramEsts(1); % (2)
                obj.Weibull_N = 1.0 / power(obj.Weibull_scale, obj.Weibull_shape);
                % Perform Kolmogorov-Smirnov test
                [obj.Weibull_h, obj.Weibull_p] = kstest(data, 'CDF', makedist('Weibull', 'a', obj.Weibull_shape, 'b', obj.Weibull_scale));
                % if h == 0           disp('The null hypothesis that the data follows a Weibull distribution cannot be rejected.');
                % else                disp('The null hypothesis that the data follows a Weibull distribution is rejected.');

                [obj.WeibullR_shape, obj.WeibullR_scale, obj.WeibullR_N, obj.WeibullR_R2, obj.WeibullR_p, obj.WeibullR_X, obj.WeibullR_Y, obj.WeibullR_xreg, obj.WeibullR_yreg] = WeibullAnalysis(obj.vals, min4Weibull);
            end        
            % normal:
            if (doGauss)
                data = obj.vals;
                pd = fitdist(data, 'Normal');
                obj.Gauss_mu = pd.mu;
                obj.Gauss_sigma = pd.sigma;
                [obj.Gauss_h, obj.Gauss_p] = kstest(data, 'CDF', makedist('Normal', 'mu', obj.Gauss_mu, 'sigma', obj.Gauss_sigma));
            end
            objout = obj;
        end
        function printHeader(obj, fid)
            fprintf(fid, 'minV\t');
            fprintf(fid, 'minLoc\t');
            fprintf(fid, 'maxV\t');
            fprintf(fid, 'maxLoc\t');
            fprintf(fid, 'meanV\t');
            fprintf(fid, 'sdivV\t');
            fprintf(fid, 'log10(sdivV)\t');
            fprintf(fid, 'covV\t');
            fprintf(fid, 'log10(covV)\t');
            fprintf(fid, 'spanV\t');
            fprintf(fid, 'nspanV\t');
    
            fprintf(fid, 'WeibullR_shape\t');
            fprintf(fid, 'Weibull_shape\t');
            fprintf(fid, 'WeibullR_scale\t');
            fprintf(fid, 'Weibull_scale\t');
            fprintf(fid, 'log10(WeibullR_scale)\t');
            fprintf(fid, 'log10(Weibull_scale)\t');
            fprintf(fid, 'WeibullR_N\t');
            fprintf(fid, 'Weibull_N\t');
            fprintf(fid, 'log10(WeibullR_N)\t');
            fprintf(fid, 'log10(Weibull_N)\t');
            
            fprintf(fid, 'WeibullR_R2\t');
            fprintf(fid, 'Weibull_h\t');
            fprintf(fid, 'Weibull_p\t');

            fprintf(fid, 'Gauss_mu\t');
            fprintf(fid, 'Gauss_sigma\t');
            fprintf(fid, 'log10(Gauss_sigma)\t');
            fprintf(fid, 'Gauss_h\t');
            fprintf(fid, 'Gauss_p\n');
        end
        function printVals(obj, fid, printHeader)
            if (nargin < 3)
                printHeader = 0;
            end
            if (printHeader)
                obj.printHeader(fid);
            end
            fprintf(fid, '%g\t', obj.minV);
            fprintf(fid, '%d\t', obj.minLoc);
            fprintf(fid, '%g\t', obj.maxV);
            fprintf(fid, '%d\t', obj.maxLoc);
            fprintf(fid, '%g\t', obj.meanV);
            fprintf(fid, '%g\t', obj.sdivV);
            fprintf(fid, '%g\t', log10(obj.sdivV));
            fprintf(fid, '%g\t', obj.covV);
            fprintf(fid, '%g\t', log10(obj.covV));
            fprintf(fid, '%g\t', obj.spanV);
            fprintf(fid, '%g\t', obj.nspanV);
    
            fprintf(fid, '%g\t', obj.WeibullR_shape);
            fprintf(fid, '%g\t', obj.Weibull_shape);
            fprintf(fid, '%g\t', obj.WeibullR_scale);
            fprintf(fid, '%g\t', obj.Weibull_scale);
            fprintf(fid, '%g\t', log10(obj.WeibullR_scale));
            fprintf(fid, '%g\t', log10(obj.Weibull_scale));
            fprintf(fid, '%g\t', obj.WeibullR_N);
            fprintf(fid, '%g\t', obj.Weibull_N);
            fprintf(fid, '%g\t', log10(obj.WeibullR_N));
            fprintf(fid, '%g\t', log10(obj.Weibull_N));
            
            fprintf(fid, '%g\t', obj.WeibullR_R2);
            fprintf(fid, '%g\t', obj.Weibull_h);
            fprintf(fid, '%g\t', obj.Weibull_p);
            
            fprintf(fid, '%g\t', obj.Gauss_mu);
            fprintf(fid, '%g\t', obj.Gauss_sigma);
            fprintf(fid, '%g\t', log10(obj.Gauss_sigma));
            fprintf(fid, '%g\t', obj.Gauss_h);
            fprintf(fid, '%g', obj.Gauss_p);
        end
    end
end