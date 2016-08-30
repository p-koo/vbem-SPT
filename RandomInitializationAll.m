function [a0 b0 mu0 sigma0 P0 D0 S2] = RandomInitializationAll(K,D,trajectory,D_cve,sigma_cve,dt)

dim = size(trajectory(1).deltaX{1},2);

% CDF of diffusivities to generate a good random initialization
bin = 1000;
MIN = 0;
MAX = max(D_cve);
edges = MIN:(MAX-MIN)/bin:MAX;
DN = histc(D_cve,edges);
cumprob = cumsum(DN)/sum(DN);

status = 1;
while status == 1

    % random population fractions greater than 5%
    P0 = ones(1,K)/K;
    for i = 1:K
        P0(i) = rand;
        while P0(i) < 0.02
            P0(i) = rand;
        end
    end
    P0 = P0/sum(P0);

    % mean location of each random population fraction along CDF 
    Dpoints = cumsum([0 P0(1:end-1)] + diff([0 P0])/2);
%     Dpoints = 1/K/2:1/K:1;
    
    % find which diffusivity corresponds to each random population fraction
    D0 = zeros(1,K);
    for i = 1:K
        index = find(cumprob > Dpoints(i),1,'first');
        D0(i) = max(edges(index),1e-6);
    end
    S2 = ones(1,K)*max(mean(sigma_cve),0);
    D0
    P0
    
    % initial mu
    mu0 = zeros(K,1);

    % initial covariance matrices
    sigma0 = zeros(D,D,K);
    for i = 1:K
        diagterm = 2*D0(i)*dt + 2*S2(i) - 2/3*D0(i)*dt;
        corrterm = 1/3*D0(i)*dt - S2(i);
        sigma0(:,:,i) = toeplitz([diagterm corrterm zeros(1,D-2)]);
    end

    a0 = InitializeTransitionMatrix(K,.97);

    b0 = EmissionProbabilitiesAll(trajectory,mu0,sigma0);
    
    if sum(sum(isnan([b0{:}]),2)) == 0 & sum(diff(D0)<1e-4) == 0
        status = 0;
    end
end




