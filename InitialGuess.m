function [Nk Sk] = InitialGuess(K,minD,maxD,sigma,random,trackparams)

dt = trackparams.dt;
R = trackparams.R;
D = trackparams.D;
N = trackparams.N;

if random == 1
    Drange = sort(rand(1,K)*maxD);
    Drange = Drange - min(Drange) + minD;
else
    Drange = minD:maxD/K:maxD;
end
sig2 = ones(1,K)*sigma^2;
disp(['Initial guess for D: ' num2str(Drange)]);

Sk = cell(1,K);
for k = 1:K
    diagterms = 2*Drange(k)*dt + 2*sig2(k) - 4*R*Drange(k)*dt;
    corrterms = 2*R*Drange(k)*dt - sig2(k);
    Sk{k} = toeplitz([diagterms corrterms  zeros(1,D-2)]);
end
Nk = ones(1,K)*round(N/K);