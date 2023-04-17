alphas = {'1.95', '1.5', '1', '0.5', '0.1'};
%alphas = {'0.1'};
%alphas = {'1.95'};
nAlpha = length(alphas);

betas = {'0.1', '0.5', '1', '1.5', '1.9'};
%betas = {'0.1', '0.5'};
nBetas = length(betas);

for ai = 1:nAlpha
    alpha = alphas{ai};
    for bi = 1:nBetas
        beta = betas{bi};
        TranslateCDFiles(alpha, beta);
%        MapCDFiles2Triangular(alpha, beta);
    end
end