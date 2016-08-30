function [Nk Sk] = SufficientStatistics(rnk,trackparams)

K = size(rnk,2);
D = trackparams.D;
vacf = trackparams.vacf;
dim = trackparams.dim;

% N_k (Eq: 10.51)
Nk = sum(rnk); 
Nk = Nk + eps;

% S_k (Eq: 10.53)
T = size(vacf,2);
diagterms = zeros(T,K);
for i = 1:T
    for j = 1:dim
        diagterms(i,:) = diagterms(i,:) + sum(rnk.*(vacf(:,i,j)*ones(1,K)))./Nk;
    end
    diagterms(i,:) = diagterms(i,:)/dim;
end

Sk = cell(1,K);
for k = 1:K
    Sk{k} = toeplitz([diagterms(:,k)' zeros(1,D-T)]);
end









