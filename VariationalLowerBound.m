function L = VariationalLowerBound(rnk,Nk,Sk,priors,est,E,trackparams)
 
D = trackparams.D;
K = size(rnk,2);

alpha0 = priors.alpha0;
invW0 = priors.invW0;
v0 = priors.v0;

alphak = est.alphak;
Wk = est.Wk;
vk = est.vk;

E_logDetLambdak = E.E_logDetLambdak;
E_logPik = E.E_logPik;
% E_xlabmdakx = E.E_xlabmdakx;

% Lower bound: E[log(P(X|Z,mu,Lambda))] (Eq. 10.71)
trSW = zeros(1,K);
for k = 1:K        
    trSW(k) = trace(Sk{k}*Wk{k});
end
E_logPX = 0.5*sum(Nk.*(E_logDetLambdak - vk.*trSW - D*log(2*pi)));%E_xlabmdakx)); %

% Lower bound: E[log(P(Z|pi))] (Eq. 10.72)
E_logPZ = sum(Nk.*E_logPik);

% Lower bound: E[log(P(X|Z,mu,Lambda))] (Eq. 10.73)
logCalpha0 = gammaln(K*alpha0) - K*gammaln(alpha0);  %(Eq. B.23)
E_logPpi = logCalpha0 + (alpha0 - 1)*sum(E_logPik);
  
% Lower bound: E[log(P(mu,Lambda))] (Eq. 10.74)        
trinvW0W = zeros(1,K); 
for k = 1:K
    trinvW0W(k) = trace(invW0*Wk{k});
end
factor = 0.5*(v0 + 1 - (1:D));
logB0 = (v0/2)*log(det(invW0)) - (v0*D/2)*log(2) ...
      - (D*(D-1)/4)*log(pi) - sum(gammaln(factor)); % (Eq. B.79)
logPLambda = K*logB0 + 0.5*(v0-D-1)*sum(E_logDetLambdak) - 0.5*sum(vk.*trinvW0W);
  
% Lower bound: E[log(Q(Z))] (Eq. 10.75)
E_logQZ = sum(sum(rnk.*log(rnk))); 

% Lower bound: E[log(P(X|Z,mu,Lambda))] (Eq. 10.71)
logCalpha = gammaln(sum(alphak)) - sum(gammaln(alphak));  %(Eq. B.23)
E_logQpi = sum((alphak-1).*E_logPik) + logCalpha;

% Lower bound: E[log(P(X|Z,mu,Lambda))] (Eq. 10.77)
H = zeros(1,K); % sum(H(q(Lamba(k)))) entropy of Wishardt distribution (Eq. B.82)
for k = 1:K
    factor = 0.5*(vk(k) + 1 - (1:D));
    logBk = (vk(k)/2)*log(det(Wk{k})) - (vk(k)*D/2)*log(2)...
            - (D*(D-1)/4)*log(pi) - sum(gammaln(factor));
    H(k) = - logBk - 0.5*(vk(k)-D-1)*E_logDetLambdak(k) + 0.5*vk(k)*D;
end
E_logQLambda = sum(.5*E_logDetLambdak - .5*D*K - H);

% Variational Lower Bound (Eq. 10.70)
L =   E_logPX + E_logPZ + E_logPpi + logPLambda - E_logQZ - E_logQpi - E_logQLambda;

%log(factorial(K)) +

