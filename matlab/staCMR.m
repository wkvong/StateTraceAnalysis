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

if ~iscell(data)
    y = outSTATS (data); % assumes general format
elseif isstruct(data{1})
    y = data; % if structured then already in stats form
else
    y = staSTATS(data); % otherwise assume within-subjects data and get stats
end

[x f Eprime] = CMRv4 (y, E);
f(find(f<tol))=0;