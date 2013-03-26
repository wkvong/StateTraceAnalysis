function [x fit exitflag output] = MR1 (y, w, E)
% Uses lsqlin to solve standard monotonic regression problem
% y is a vector of values - the data
% MR conducts monotonic regression on each column separately
% w is the weight matrix for the (lsqlin) fit: either a vector of positive
% numbers or a positive definite matrix
% E is an adjacency matrix coding the partial order model
% 
% returns:
% x = best fitting values
% Ffit = weighted least squares fit
% exitflag = vector of exit flags for fits for each variable

n = numel(y);
if size(y,1)==1 % convert to column vector
    y=y';
end
if nargin==1
    w = eye(n);
    E = zeros(n,n);
elseif nargin==2
    E = zeros(n,n);
end
if isempty(w)
    w = eye(n);
end
if isempty(E)
    E = zeros(n,n);
end
if iscell(E)
    adj = cell2adj(1:numel(y),E);
else
    adj = E;
end
if sum(sum(adj)) > 0
    A = adj2ineq (adj); % John's code to turn an adjacency model to a set of inequalities
    b = zeros(size(A,1),1);
else
    A = []; b =[];
end

options = optimset ('LargeScale','off', 'display','off'); % turn off display and use medium algorithm to avoid warning
% do monotonic regression for each column of y
if isvector(w)
    C = diag(sqrt(w));
else
    C = sqrtm(w);
end

d = C*y;

x0 = repmat(mean(y),n,1); % starting point is a vector of means

[x fit resid exitflag output] = lsqlin(C,d,A,b,[],[],[],[],x0,options);
%fit = real(fit)/(numel(y)-1);