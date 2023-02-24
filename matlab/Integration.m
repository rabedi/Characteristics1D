function integrationVals = Integration(times, integrandVals)
numTimes = length(times);
delTime = times(2:numTimes) - times(1:numTimes - 1);
% fact(1) = 0.5 * delTime(1);
% for i = 2:numTimes - 1
%     fact(i) = 0.5 * (delTime(i - 1) + delTime(i));
% end
% fact(numTimes) = 0.5 * delTime(numTimes - 1);
% 
% integrationVals = zeros(m, n);
[m, n] = size(integrandVals);
integrationVals = zeros(m, n);
for j = 1:n
    integrationVals(1, j) = 0.0;
end
for i = 2:numTimes
    hdelT = 0.5 * delTime(i - 1);
    for j = 1:n
        integrationVals(i, j) = integrationVals(i - 1, j) + hdelT * (integrandVals(i - 1, j) + integrandVals(i, j));
    end
end