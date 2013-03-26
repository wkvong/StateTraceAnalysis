function [x f Eprime] = staCMR (data, E)
% data is cell array of data or structured output from staSTATS
% E is partial order
% returns:
% x = best fitting CMR values to y-means
% f = fit statistic
tol=10e-5;
if nargin == 1
    E = {};
end
if isempty(E)
    E = {};
end
y = data;
if ~iscell(y)
    y = {y};
end
if ~isstruct(y{1})
    y = staSTATS(y);
end
[x f Eprime] = CMRv4 (y, E);
f(find(f<tol))=0;
if ~iscell(data)
    x = x{1};
    f = f(1);
end
