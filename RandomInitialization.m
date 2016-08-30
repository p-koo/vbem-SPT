function [D0 S2]= RandomInitialization(K,deltaX,D_cve,sigma_cve)

% bin size dimensions
D = size(deltaX{1},1);

% CDF of diffusivities to generate a good random initialization
bin = 1000;
MIN = 0;
MAX = max(D_cve);
edges = MIN:(MAX-MIN)/bin:MAX;
DN = histc(D_cve,edges);
cumprob = cumsum(DN)/sum(DN);

% random population fractions greater than 2%
P0 = ones(1,K)/K;
for i = 1:K
    P0(i) = rand;
    while P0(i) < 0.02
        P0(i) = rand;
    end
end
P0 = P0/sum(P0);

% mean location of each random population fraction along CDF 
Dpoints = cumsum([0 P0(1:end-1)] + diff([0 P0])*.5);

% find which diffusivity corresponds to each random population fraction
D0 = zeros(1,K);
for i = 1:K
    index = find(cumprob > Dpoints(i),1,'first');
    D0(i) = max(edges(index),1e-6);
end
disp(['D0 = ' num2str(D0)]);

% initial square localization noise 
S2 = ones(1,K)*max(mean(sigma_cve),0);


