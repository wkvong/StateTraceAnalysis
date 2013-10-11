function [p datafit fits] = CMRMNfits (nsample, data, E)
% Parametric bootstrap sampling test of multinomial data
% data is NCOND x NRESP matrix of counts
% E is partial order in form of cell array or adjacency matrix

if nargin==2
    E={};
end
y = data;
[x f2] = staMRMN (y, E);
[x f1] = staCMRMN (y, E);
datafit = [f1 f2 f1-f2];

rng ('shuffle');
fits = zeros(nsample,3);
parfor isample=1:nsample
    % bootstrap sample
    yb = resampleMN (y); 
    % fit model to bootstrap data
    x = staCMRMN (yb, E);
    % resample model
    yr = resampleMN (x); 
    % fit model
    [x f2] = staMRMN (yr, E);
    [x f1] = staCMRMN (yr, E);
    fits(isample,:) = [f1 f2 f1-f2];
end

% calculate p
k = find(fits(:,3) >= datafit(3));
p = numel(k)/nsample;


function yr = resampleMN(x)
s = sum(x,2); n = round(s);
p = x./repmat(s,1,size(x,2));
yr = mnrnd (n,p);

