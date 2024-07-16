function [X, Y, xreg, yreg, R2, shape, scale, N] = loadWeibullPlotData(fileNWOExt)
fn = [fileNWOExt, '.mat'];
% Load the structure from the binary file
loaded_data = load(fn);

% Access the variables from the loaded structure
WeibullR_struct_loaded = loaded_data.WeibullR_struct;

% Retrieve individual variables
p = WeibullR_struct_loaded.p;
xreg = WeibullR_struct_loaded.xreg;
yreg = WeibullR_struct_loaded.yreg;
X = WeibullR_struct_loaded.X;
Y = WeibullR_struct_loaded.Y;
R2 = WeibullR_struct_loaded.R2;
shape = WeibullR_struct_loaded.shape;
scale = WeibullR_struct_loaded.scale;
N = WeibullR_struct_loaded.N;



