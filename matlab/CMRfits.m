function [p datafit fits] = CMRfits (nsample, data, E, E1)
% if E1 is specified then compares non-overlap model with 1D model
% else compares 1D model to 2D model
% computes empirical p-value for hypothesis that means across NCOND
% conditions are monotonically increasing
% nsample = no. of Monte Carlo samples (about 10000 is good)
% data = NVAR cell array of NSUBJ x NCOND matrices of observations
% E = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
% condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
% returns:
% p = empirical p-value
% datafit = observed fit of monotonic (1D) model
% fits = nsample vector of fits of Monte Carlo samples (it is against this
% distribution that datafit is compared to calculate p)
% Note: last column of outputs correspond to sum of fits
%
if ~iscell(data)
    type = 1; % data in "general format" i.e. one line per subject
    nvar = numel(unique(data(:,3)));
else
    type = 0; % data in within-subject cell array format
    nvar = numel(data);
end
if nargin < 3
    E = {};
    E1 = {};
elseif nargin == 3
    E1 = {};
end

if type == 1;
    ys = outSTATS(data);
else
    ys = staSTATS(data); 
end

if isempty(E1)
    if ~isempty(E)
        [x2 f2] = staMR (ys, E);
    else
        f2 = 0;
    end
    [x1 f1] = staCMR (ys, E);
else
    [x2 f2] = staCMR (ys, E);
    [x1 f1] = staCMR (ys, E1);
end
f = f1 - f2; datafit = [f sum(f)];

%rng ('shuffle');
rand('seed',sum(100*clock));
fits = zeros(nsample,nvar);
parfor isample=1:nsample
    % bootstrap sample
    yb = bootstrap (data, type);
    % fit 1D model to bootstrap data
    if type==0
        y = staSTATS(yb);
    else
        y = outSTATS(yb);
    end
    if isempty(E1)
        [x f] = staCMR (y, E);
    else
        [x f] = staCMR (y, E1);
    end
    % resample model
    yr = resample (x, y, type);
    % fit 1D and 2D models
    %y = staSTATS(yr);
    y = yr;
    if isempty(E1)
        if ~isempty(E)
            [x2 f2] = staMR (y, E);
        else
            f2 = 0;
        end
        [x1 f1] = staCMR (y, E);
    else
        [x2 f2] = staCMR (y, E);
        [x1 f1] = staCMR (y, E1);
    end
    f = f1 - f2; 
    fits(isample,:) = f; % store Monte Carlo fits
end
fits = [fits sum(fits,2)];  % add sum column
% calculate p
p = zeros(size(datafit));
for i=1:size(fits,2)
    k = find(fits(:,i) >= datafit(i));
    p(i)=numel(k)/nsample;
end

function yb = bootstrap (y, type)
% draws bootstrap sample from data
if type==0 % y = cell array of NSUBJ x NCOND matrices
    for ivar = 1:numel(y)
        a = y{ivar};
        nsub = size(a,1);
        yy = zeros(size(a)); b = zeros(1,nsub);
        v=nsub;
        while v==nsub
            for isub = 1:nsub
                u = randperm(nsub);
                yy(isub,:) = a(u(1),:);
                b(isub)=u(1);
            end
            v = sum(b==b(1)); % check if only sampled a single subject
        end
        yb{ivar}=yy;
    end
else % y = "general format"
    cond = unique (y(:,2));
    var = unique (y(:,3));
    yb = y;
    for j=1:numel(var)
        for i=1:numel(cond)
            k = find(y(:,2)==cond(i) & y(:,3)==var(j));
            a = y(k,:);
            r = floor(rand(1,numel(k))*numel(k))+1;
            yb(k,:) = a(r,:);
        end
    end
end

function yr = resample(x, y, type)
yr = y;
for ivar = 1:numel(x)
    if type==0
        sigma = y{ivar}.cov./y{ivar}.n;
    else
        sigma = zeros(size(y{ivar}.cov));
        k = find(y{ivar}.n > 0);
        sigma(k)=y{ivar}.cov(k)./y{ivar}.n(k);
    end
    yr{ivar}.means = mvnrnd (x{ivar}, sigma)';
end