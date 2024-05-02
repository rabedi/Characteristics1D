function plot_epsDot_conv_plots()
% AT1 and 2 use pPrimes so that sigmaMax = 1
CZM_normalization4AT = 1;
% bPrime of AT1 and AT2 are 1
CZM_normalization4AT = 0;


plot3 = 0; % plot all combinations [AT1, EPH, exHR] in one plot
plot2 = 0; % only two things in one plot
plot1 = 0; % only one thing changes in each plot
plot_df_comp = 1;

gvid = genVisIndexDat;
gvid.colNamesLatex{1} = '$$ \mathrm{log}({\bar{\sigma}}_M) $$';
gvid.colNamesLatex{2} = '$$ \mathrm{log}({\bar{\epsilon}}_M) $$';
gvid.colNamesLatex{3} = '$$ \mathrm{log}({\bar{\epsilon}}_f) $$';
gvid.colNamesLatex{4} = '$$ \mathrm{log}({\bar{\phi}}_{f}) $$';
gvid.colNamesLatex{5} = '$$ \mathrm{log}({\bar{\phi}}_{u}) $$';
gvid.colNamesLatex{6} = '$$ \mathrm{log}({\bar{\phi}}_{M}) $$';
gvid.colNamesLatex{7} = '$$ B^\phi $$';
gvid.colNamesLatex{8} = '$$ B^\epsilon $$';

gvid.xLabelLatex = '$$ \mathrm{log}({\bar{\dot{\epsilon}}}_0) $$';
% group 1: reserve -> not used
gvid.groupNames{1} = 'none';
% group 2
df_ss = {'0.01', '0.1', '1', '10', '100'};
n_df_ss = length(df_ss);
gvid.groupNames{2} = df_ss;
% group 3: 0D vs 0D HR
eHR = {'0D', '0D (HRA)'};
eHRns = {'0D_ex', '0D_HRA'};
gvid.groupNames{3} = eHR;
% group 4: EPF, PPF, HPF
EPH = {'EPF', 'PPF', 'HPF'};
gvid.groupNames{4} = EPH;
% group 5: model name
model_ss = {'AT1', 'AT2', 'CZM-W'};
num_model = length(model_ss);
gvid.groupNames{5} =  model_ss;


offset = 0;
offsetConv = 0;
if (CZM_normalization4AT == 0)
    offset = 5;
    offsetConv = 500;
end

PHnames = {'P', 'H'};
i1 = 1;
i2 = 1;



for i2 = 1:n_df_ss
    df_s = df_ss{i2};
    df_v(i2) = str2num(df_s);
    for option = 31:32
        option_i = option - 30;
        for isHyper = 0:1
            hi = isHyper + 1;
            for mi = 1:num_model
                model_s = model_ss{mi};
                forlderName = ['option_', num2str(option + offset)];
%                forlderName = [forlderName, '_', PHnames{hi}, '_', model_s, '_lcRat_none_df_none'];
                forlderName = [forlderName, '_', PHnames{hi}, '_', model_s, '_lcRat_none_df_', df_s];
                fileName = [forlderName, '/scalars.csv'];
                [gvid.data{i1}{i2}{option_i}{hi + 1}{mi}, gvid.xVals{i1}{i2}{option_i}{hi + 1}{mi}] = Read_PF_cvs(fileName);
            end
        end
    end
end
ind_df1 = 3;
added_indices = cell(0);
added_indices_mi = cell(3);
hi = 1; % EPF
for mi = 1:num_model
    model_s = model_ss{mi};
    fileName = ['model-', model_s];
    [gvid.data{1}{ind_df1}{1}{hi}{mi}, gvid.xVals{1}{ind_df1}{1}{hi}{mi}] = Read_PF_cvs(fileName);
    indv = [1, ind_df1, 1, hi, mi];
    added_indices{mi} = indv;
    added_indices_mi{mi}{1} = indv;
end

fclose('all');


starts = [1, ind_df1, 1, 2, 1];
ends = [1, ind_df1, 2, 3, num_model];
addedFolderFileName = '';
% plotting everything
optionNo = 110 + offsetConv;
clrPos = 5;
lineStylePos = 3;
markerPos = 4;
if (plot3)
    gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices)
end

if (plot2) 
    %%%%%%%%%%%%%%
    % 1 -> cDRat, 2 -> df, 3 -> [1 exact, 2 HR], 4 -> EPH, 5 -> AT1, ...
    % 3 ways of combining 2 things in one plot
    optionNo = 120 + offsetConv;
    % EPH + exact, HR in one plot
    clrPos = 4;
    lineStylePos = 3;
    markerPos = -1;
    for mi = 1:num_model
        model_s = model_ss{mi};
        addedFolderFileName = ['model', num2str(mi), '_', model_s];
        starts = [1, ind_df1, 1, 2, mi];
        ends = [1, ind_df1, 2, 3, mi];
        gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices_mi{mi});
    end
    
    optionNo = 121 + offsetConv;
    % model(AT,...) + exact, HR in one plot
    clrPos = 5;
    lineStylePos = 3;
    markerPos = -1;
    for hi = 2:3
        eph_s = EPH{hi};
        addedFolderFileName = ['EPH', num2str(hi), '_', eph_s];
        starts = [1, ind_df1, 1, hi, 1];
        ends = [1, ind_df1, 2, hi, num_model];
        gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices);
    end


    optionNo = 122 + offsetConv;
    % EPH + model(AT,...) 
    clrPos = 5;
    lineStylePos = 4;
    markerPos = -1;
    for option_i = 1:2
        option_s = eHRns{option_i};    
        starts = [option_i, ind_df1, 1, 2, 1];
        ends = [option_i, ind_df1, 2, 3, num_model];
        addedFolderFileName = ['ex_hr', num2str(option_i), '_', option_s];
        starts = [1, ind_df1, option_i, 2, 1];
        ends = [1, ind_df1, option_i, 3, num_model];
        gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices);
    end
end


if (plot1) 
    %%%%%%%%%%%%%%
    % 1 -> cDRat, 2 -> df, 3 -> [1 exact, 2 HR], 4 -> EPH, 5 -> AT1, ...
    % 3 ways of having only one thing in a plot
   optionNo = 130 + offsetConv;
    % model only changed on plot
    clrPos = 5;
    lineStylePos = -1;
    markerPos = -1;
    for hi = 2:3
        eph_s = EPH{hi};
        addedFolderFileName1 = ['EPH', num2str(hi), '_', eph_s];

        for option_i = 1:2
            option_s = eHRns{option_i};    
            addedFolderFileName2 = ['ex_hr', num2str(option_i), '_', option_s];

            addedFolderFileName = [addedFolderFileName1, '_', addedFolderFileName2];
            starts = [1, ind_df1, option_i, hi, 1];
            ends = [1, ind_df1, option_i, hi, num_model];
%            gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices);
        end
    end


   optionNo = 131 + offsetConv;
    % 1 -> cDRat, 2 -> df, 3 -> [1 exact, 2 HR], 4 -> EPH, 5 -> AT1, ...
    % EPH only changed in plot
    clrPos = 4;
    lineStylePos = -1;
    markerPos = -1;

    for mi = 1:num_model
        model_s = model_ss{mi};
        addedFolderFileName1 = ['model', num2str(mi), '_', model_s];

        for option_i = 1:2
            option_s = eHRns{option_i};    
            addedFolderFileName2 = ['ex_hr', num2str(option_i), '_', option_s];

            addedFolderFileName = [addedFolderFileName1, '_', addedFolderFileName2];
            starts = [1, ind_df1, option_i, 2, mi];
            ends = [1, ind_df1, option_i, 3, mi];
            gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices_mi{mi});
        end
    end
end

if (plot_df_comp == 0)
    return;
end


optionNo = 200 + offsetConv;
clrPos = 2;
lineStylePos = -1;
markerPos = -1;
for option_i = 1:2
    option_s = eHRns{option_i};    
    addedFolderFileName1 = ['x', num2str(option_i), '_', option_s];
    for mi = 1:num_model
        model_s = model_ss{mi};
        addedFolderFileName2 = ['m', num2str(mi), '_', model_s];
        for hi = 2:3
            eph_s = EPH{hi};
            addedFolderFileName3 = ['p', num2str(hi), '_', eph_s];
            addedFolderFileName = [addedFolderFileName1, '_', addedFolderFileName2, '_', addedFolderFileName3];
            starts = [1, 1, option_i, hi, mi];
            ends = [1, n_df_ss, option_i, hi, mi];
            gvid.Visualize(optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices_mi{mi});
        end
    end
end