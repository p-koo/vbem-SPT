function  L = LogMultiVariateNormal2(deltax,mu,C)
 

D = length(deltax);

T = cholcov(C);
if ~isempty(T)
    deltax2 = inv(T)'*(deltax-mu);
    logdetC = 2*log(det(T));
    L = -D/2*log(2*pi) -.5*logdetC - .5*deltax2'*deltax2;
else
    threshold = 1e-8;
    [vec val] = eig(C);
    val = diag(val);
    val(val < threshold) = threshold;
    logdetC = sum(log(val));
    invL = diag(1./val);
    deltax2 = vec'*(deltax-mu);
    L = -D/2*log(2*pi) -.5*logdetC - .5*deltax2'*invL*deltax2;
end

