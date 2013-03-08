function [x f exitflag] = staMR (data, E)
% fits monotonic regression model to data according to partial order
% data is cell array of data or structured output from staSTATS
% E is partial order
% returns:
% x = best fitting MR values to y-means
% f = fit statistic
tol = 10e-5;
if nargin == 1
    E = {};
end
if ~iscell(E) && isvector(E) % convert to cell array if only lazy vector is specified
    E = {E};
end
y = data;
if ~iscell(y)
    y = {y};
end
if ~isstruct(y{1})
    y = staSTATS(y);
end
x = {}; f=zeros(1,numel(y)); exitflag=zeros(size(f));
for i = 1:numel(y)
    [x{i} f(i) exitflag(i)] = MR1 (y{i}.means, y{i}.weights, E);
end
f(find(f<tol))=0;

if ~iscell(data)
    x = x{1};
end

