function E = UsefulExpectations(estimates,deltax,trackparams)

D = trackparams.D;
N = trackparams.N;
dim = trackparams.dim;
C = trackparams.C;

alphak = estimates.alphak;
Wk = estimates.Wk;
vk = estimates.vk;
K = length(alphak);

% E-step: xlambdakx (Eq. 10.64)
E_xlabmdakx = zeros(N,K);
for k = 1:K
    for n = 1:N
%         E_xlabmdakx(n,k) =  vk(k)*trace(Wk{k}*C{n});
        for j = 1:dim
            E_xlabmdakx(n,k) =  E_xlabmdakx(n,k) + vk(k)*(deltax{n}(:,j)'*Wk{k}*deltax{n}(:,j));
        end
        E_xlabmdakx(n,k) = E_xlabmdakx(n,k)/dim;
    end
end

% E-step: logLambdaTilde (Eq. 10.65)
E_logDetLambdak = zeros(1,K);
for k = 1:K    
    factor = (vk(k)+1-(1:D))/2;
    sumPsi = sum(psi(0,factor));    
    logdetWk = LogDeterminant(Wk{k});
    E_logDetLambdak(k) = sumPsi + D*log(2) + logdetWk;
end

% E-step: logPiTilde (Eq. 10.66)
E_logPik = psi(0,alphak) - psi(0,sum(alphak));

E.E_xlabmdakx = E_xlabmdakx;
E.E_logDetLambdak = E_logDetLambdak;
E.E_logPik = E_logPik;

