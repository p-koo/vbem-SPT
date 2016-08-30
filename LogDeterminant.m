function [logdetC invC]= LogDeterminant(S)

threshold = 1e-8;
% try
    [vec val] = eig(S);
    val = diag(val);
    index = find(val <=threshold);
    val(index) = threshold;
    logdetC = sum(log(val));
    invC = vec*diag(1./val)*vec';
% catch
%     [L U P] = lu(S);
%     du = diag(U);
%     index = find(du <=threshold);
%     du(index) = threshold;
%     c = abs(det(P));
%     logdetC = log(c) + sum(log(abs(du)));
%     invC = inv(U)*inv(L);

end