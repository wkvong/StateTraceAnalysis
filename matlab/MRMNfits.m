function [p datafit fits] = MRMNfits (nsample, data, E)
% Parametric bootstrap sampling test of multinomial data
% data is NCOND x NRESP matrix of counts
% E is partial order in form of cell array or adjacency matrix

if nargin==2
    E={};
end
y = data;
[x f] = staMRMN (y, E);
datafit = f;

rng ('shuffle');
fits = zeros(nsample,1);
parfor isample=1:nsample
    % bootstrap sample
    yb = bootstrap (y);
    % fit model to bootstrap data
    x = staMRMN (yb, E);
    % resample model
    yr = resample (x, yb);
    % fit model
    [x f] = staMRMN (yr, E);
    fits(isample) = f;
end

% calculate p
k = find(fits >= datafit);
p = numel(k)/nsample;

function yr = resample(x,y)
n = sum(y,2);
p = x./repmat(n,1,size(x,2));
p(:,size(p,2))=1-sum(p(:,1:size(p,2)-1),2); % round off sum to 1
yr = mnrnd (n,p);


    