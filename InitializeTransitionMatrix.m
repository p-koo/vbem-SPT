function a = InitializeTransitionMatrix(K,A)

% initialize transmission matrix
a = diag(ones(1,K)*A) + ones(K,K)*(1-A);
a = a./(sum(a,2)*ones(1,K));
