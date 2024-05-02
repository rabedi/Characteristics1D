option = 30;
isHyper = -1;
isHyper = 0;
isApproximate = 2;
% AT1, AT2, CZM-W, CZM-L are the options
model_s = 'AT1';
la_s = 'none';
l_cD2c_s = 'none';
df_s = 'none';
CZM_normalization4AT = -1;
bTimesPiCZM_s = '1.0';
CZM_modelName = 'Linear';

% change to 0 for bPrime = 1 for AT1 and AT2
CZM_normalization4AT = 0;


model_ss = {'AT1', 'AT2', 'CZM-W'};

if 0
    df_ss = {'0.01', '0.1', '1', '10', '100'};
    sz_df = length(df_ss);
    for dfi = 1:sz_df
        df_s = df_ss{dfi};
        df = str2num(df_s);
        for option = 30:32
            for isHyper = 0:1
                for mi = 1:length(model_ss)
                    model_s = model_ss{mi};
                    fprintf(1, 'df = %s, option = %d, hyper = %d, model, %s\n', df_s, option, isHyper, model_s);
                    plot_print_PFs(option, isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName);
                end
            end
        end
    end
end


if 0
    for option = 30:32
        for isHyper = 0:1
            for mi = 1:length(model_ss)
                model_s = model_ss{mi};
                fprintf(1, 'option = %d, hyper = %d, model, %s\n', option, isHyper, model_s);
                plot_print_PFs(option, isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName);
            end
        end
    end
end

if 1
    options = [40, 50];
    n_options = length(options);
    la_ss = {'2', '1', '0', '-1'};
    n_la_ss = length(la_ss);
    for li = 1:n_la_ss
        la_s = la_ss{li};
        for oi = 1:n_options
            option = options(oi);
            for isHyper = 0:1
                for mi = 1:length(model_ss)
                    model_s = model_ss{mi};
                    fprintf(1, 'la = %s, option = %d, hyper = %d, model, %s\n', la_s, option, isHyper, model_s);
                    plot_print_PFs(option, isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName);
                end
            end
        end
    end
end