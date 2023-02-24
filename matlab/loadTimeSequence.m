function [dataOut, delu, delv, sigmaL, sigmaR, diss, del_e, exists] = loadTimeSequence(fileName, dim, ringProb)
% data is indexed as follows
%time	sigma	vL		vR		uL	uR	damage	damageSource	v_r_final	v_r_source_final	sigma_theta_source_final
if (nargin < 1)
    fileName = 'xPos_1_InterfaceRawFinalSln.txt';
end
if (nargin < 2)
    dim = 1;
end
if (nargin < 3)
    ringProb = 0;
end
fid = fopen(fileName, 'r');
if (fid < 0)
    dataOut = cell(0);
    exists = 0;
    return;
end
exists = 1;
dat = load(fileName, '-ascii');
posField = 1;
pos_scalar = 1;
% time
times = dat(:, pos_scalar);
dataOut{posField} = times;
% sigma
posField = posField + 1;
pos_scalar = pos_scalar + 1;
sigmaL = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = sigmaL;
pos_scalar = pos_scalar + dim - 1;

posField = posField + 1;
pos_scalar = pos_scalar + 1;
sigmaR = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = sigmaR;
pos_scalar = pos_scalar + dim - 1;

% vL
posField = posField + 1;
pos_scalar = pos_scalar + 1;
vL = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = vL;
pos_scalar = pos_scalar + dim - 1;

% vR
posField = posField + 1;
pos_scalar = pos_scalar + 1;
vR = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = vR;
pos_scalar = pos_scalar + dim - 1;

% uL
posField = posField + 1;
pos_scalar = pos_scalar + 1;
uL = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = uL;
pos_scalar = pos_scalar + dim - 1;

% uR
posField = posField + 1;
pos_scalar = pos_scalar + 1;
uR = dat(:, pos_scalar:pos_scalar + dim - 1);
dataOut{posField} = uR;
pos_scalar = pos_scalar + dim - 1;

delu = uR - uL;

pL = sigmaL .* vL;
pR = sigmaR .* vR;
delv = vR - vL;

eL = Integration(times, pL);
eR = Integration(times, pR);
del_e = eR - eL;
[m, n] = size(del_e);
diss = zeros(m, 1);
for i = 1:m
    for j = 1:n
        diss(i) = diss(i) + del_e(i, j);
    end
end
dissBK = diss;
diss(1) = 0;
for i = 2:m
    diss(i) = diss(i - 1);
    for j = 1:n
        diss(i) = diss(i) + 0.5 * (delu(i) - delu(i - 1)) * (sigmaR(i) + sigmaR(i - 1));
    end
end
%diff = (diss - dissBK)


% delu
posField = posField + 1;
dataOut{posField} = delu;

% pL
posField = posField + 1;
dataOut{posField} = pL;

% pR
posField = posField + 1;
dataOut{posField} = pR;

% eL
posField = posField + 1;
dataOut{posField} = eL;

% del_e
posField = posField + 1;
dataOut{posField} = del_e;

% diss
posField = posField + 1;
dataOut{posField} = diss;
%ind = find(diss < -0.003)
%if (length(ind) > 0)
%    ind
%end

% damage
posField = posField + 1;
pos_scalar = pos_scalar + 1;
dataOut{posField} = dat(:, pos_scalar);

% damageSource
posField = posField + 1;
pos_scalar = pos_scalar + 1;
dataOut{posField} = dat(:, pos_scalar);

% deluEff
posField = posField + 1;
pos_scalar = pos_scalar + 1;
dataOut{posField} = dat(:, pos_scalar);

%%%%
%%%%

if (ringProb)
    % v_r_final
    posField = posField + 1;
    pos_scalar = pos_scalar + 1;
    dataOut{posField} = dat(:, pos_scalar);
    
    % v_r_source_final
    posField = posField + 1;
    pos_scalar = pos_scalar + 1;
    dataOut{posField} = dat(:, pos_scalar);
    
    % sigma_theta_source_final
    posField = posField + 1;
    pos_scalar = pos_scalar + 1;
    dataOut{posField} = dat(:, pos_scalar);
end


