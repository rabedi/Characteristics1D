function plotInitialMesh(serialNum, versionNo)
if nargin < 1
    serialNum = 0;
end
if nargin < 2
    versionNo = 4;
end

srstr = num2str(serialNum);
vstr = '';
if (versionNo >= 0)
    vstr = ['V_', num2str(versionNo), '_'];
end
root = ['../run_', vstr, srstr, '/'];
fni = [root, 'field0_initialResolution.txt'];
fnf = [root, 'field0_finalResolution.txt'];

di = load(fni);
df = load(fnf);
[m, n] = size(df)
x = -0.5:1/(m - 1):0.5;
plot(x, df);
hold on;
[m, n] = size(di)
x = -0.5:1/(m - 1):0.5;
plot(x, di);
legend({'final', 'initial'});
