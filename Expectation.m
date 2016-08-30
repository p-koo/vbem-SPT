function rnk = Expectation(E,trackparams)

N = trackparams.N;
D = trackparams.D;
K = length(E.E_logPik);

% E-step: r (Eq. 10.46)
logRho = ones(N,1)*(E.E_logPik+0.5*E.E_logDetLambdak) ...
                                     - 0.5*E.E_xlabmdakx - D/2*log(2*pi);

% E-step: rnk (Eq. 10.49)
logSumRho = logsumexp(logRho,2);
logr = logRho - repmat(logSumRho, 1,K);
rnk = exp(logr);
rnk(find(rnk(:) < eps)) = eps; %1e-100;




