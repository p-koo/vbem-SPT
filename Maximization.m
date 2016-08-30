function estimates = Maximization(Nk,Sk,priors)

% priors
alpha0 = priors.alpha0;
invW0 = priors.invW0;
v0 = priors.v0;

K = length(Sk);

D = length(Sk{1});
T = 2;

% alpha_k (Eq: 10.58)
alphak = alpha0 + Nk;

% W_k (Eq: 10.62)
invWk = cell(1,K); Wk = cell(1,K);
for k = 1:K
    invWk{k} = invW0 + Nk(k)*Sk{k};
    Wk{k} = inv(invWk{k});
%     diagterm = CalculateVACF(Wk{k},D);
%     Wk{k} = toeplitz([diagterm(1:T) zeros(1,D-T)]); 
end

% sufficent statistics: v_k (Eq: 10.63)
vk = v0 + Nk + 1;

estimates.alphak = alphak;
estimates.Wk = Wk;
estimates.vk = vk;

