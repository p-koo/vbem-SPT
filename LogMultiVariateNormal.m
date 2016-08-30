function  L = LogMultiVariateNormal(deltax,deltay,mu,C)
 
threshold = 1e-8;
N = length(deltax);

L1 = Inf; L2 = Inf;
try
    [vec val] = eig(C);
    val = diag(val);
    val(val < threshold) = threshold;
    logdetC = sum(log(val));
    invL = diag(1./val);
    
    deltax2 = vec'*(deltax-mu(1));
    L1 = -N/2*log(2*pi) -.5*logdetC - .5*deltax2'*invL*deltax2;
    deltax2 = vec'*(deltay-mu(2));
    L2 = -N/2*log(2*pi) -.5*logdetC - .5*deltax2'*invL*deltax2;    
end

L = L1 + L2;
