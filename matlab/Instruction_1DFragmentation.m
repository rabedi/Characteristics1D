classdef Instruction_1DFragmentation
    properties
        dd2 = 0.9;
        ldelc = -1.5;
        la = 1.0;
        llc = -3.0;
        shape = 2.0;
        % mesh is refined with factor resolutionFactor
        resolutionFactor = 1;
        ssoFS = 'min'; %mean_arithmetic, mean_harmonic, valStart
        % dt_resap = 1, time step remains the same as the finest mesh, 0 ->
        % time step changes based on the current mesh size
        dt_resap = 1;

        % indices
        ind_dd2 = 0;
        ind_ldelc = 0;
        ind_la = 0;
        ind_llc = 0;
        ind_shape = 0;
        ind_resolutionFactor = 0;
        ind_ssoFS = 0;
        ind_dt_resap = 0;

        success = 0;
        version = -1;
    end
    methods
        function objout = read(obj, fid, versionOffset)
            if nargin < 3
                versionOffset = 0;
            end
            objout = obj;
            objout.resolutionFactor = 1;
            objout.ind_resolutionFactor = 0;
            objout.ssoFS = 'min';
            objout.ind_ssoFS = 0;
            objout.dt_resap = 1;
            objout.ind_dt_resap = 0;
            
            buf = fscanf(fid, '%s', 1);
            status = feof(fid);
            if (status == 1)
                objout.success = 0;
                return;
            end
            objout.version = fscanf(fid, '%d', 1) + versionOffset;
            buf = fscanf(fid, '%s', 3);
            buf = fscanf(fid, '%s', 1);
            while (strcmp(buf, '}') == 0)
                name = fscanf(fid, '%s', 1);
                buf = fscanf(fid, '%s', 1);
                ind = fscanf(fid, '%d', 1);
                buf = fscanf(fid, '%s', 1);
                str = fscanf(fid, '%s', 1);
                num = str2num(str);
                if (strcmp(name, 'dd2') == 1)
                    objout.dd2 = num;
                    objout.ind_dd2 = ind;
                elseif (strcmp(name, 'ldelc') == 1)
                    objout.ldelc = num;
                    objout.ind_ldelc = ind;
                elseif (strcmp(name, 'la') == 1)
                    objout.la = num;
                    objout.ind_la = ind;
                elseif (strcmp(name, 'llc') == 1)
                    objout.llc = num;
                    objout.ind_llc = ind;
                elseif (strcmp(name, 'shape') == 1)
                    objout.shape = num;
                    objout.ind_shape = ind;
                elseif (strcmp(name, 'resolutionFactor') == 1)
                    objout.resolutionFactor = num;
                    objout.ind_resolutionFactor = ind;
                elseif (strcmp(name, 'ssoFS') == 1)
                    objout.ssoFS = str;
                    objout.ind_ssoFS = ind;
                elseif (strcmp(name, 'dt_resap') == 1)
                    objout.dt_resap = num;
                    objout.ind_dt_resap = ind;
                end
                buf = fscanf(fid, '%s', 1);
            end
            objout.success = 1;
            if (objout.success == 1)
                objout.Print();
            end
        end
        function Print(obj)
            fn = ['Instructions_V_', num2str(obj.version), '.txt'];
            fido = fopen(fn, 'w');
            fprintf(fido, '8');
            fprintf(fido, '\tdd2\t%d\t%g', obj.ind_dd2, obj.dd2);
            fprintf(fido, '\tldelc\t%d\t%g', obj.ind_ldelc, obj.ldelc);
            fprintf(fido, '\tla\t%d\t%g', obj.ind_la, obj.la);
            fprintf(fido, '\tllc\t%d\t%g', obj.ind_llc, obj.llc);
            fprintf(fido, '\tshape\t%d\t%g', obj.ind_shape, obj.shape);
            fprintf(fido, '\tresolutionFactor\t%d\t%g', obj.ind_resolutionFactor, obj.resolutionFactor);
            fprintf(fido, '\tssoFS\t%d\t%s', obj.ind_ssoFS, obj.ssoFS);
            fprintf(fido, '\tdt_resap\t%d\t%g', obj.ind_dt_resap, obj.dt_resap);
            fclose(fido);
        end
  end
end