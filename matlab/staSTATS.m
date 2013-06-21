function output = staSTATS (data)
% returns statistics on NSUB x NCOND matrix of data
% data may be a cell array of matrices
% returns means, cov, nsub and weights
% output.means = observed means
% output.n = number of subjects
% output.cov = observed covariance matrix
% output.covs = pooled covariance matriex (compound symmetry) - not used 
% output.weights = weight matrix for monotonic regression

y = data;
if ~iscell(data)
    y = {y};
end
output = cell(1,numel(y));
for idata = 1:numel(y)
    yy = y{idata};
    u = sum(isnan(yy),2); k=find(u==0); yy=yy(k,:); % delete nans
    out.means = mean(yy,1);
    out.n = size(yy,1);
    if out.n > 1
        out.cov = cov(yy);
        out.weights = out.n./diag(out.cov);
    else
        out.cov = zeros(size(yy,2),size(yy,2));
        out.weights = ones(size(yy,2),1);
    end
    out.weights=out.weights';
    out.lm = loftusmasson(yy);
 %   out.lm = diag(diag(out.cov));
    output{idata}=out;
end
if ~iscell(data)
    output = output{1};
end

function r = loftusmasson(yy)
% returns Loftus & Masson within subjects error based on Mean Square residual
y_cond = repmat(mean(yy,1), size(yy,1), 1);
y_subj = repmat(mean(yy,2), 1, size(yy,2));
y_mean = repmat(mean(mean(yy,2)), size(yy,1), size(yy,2));
ya = yy - y_cond - y_subj + y_mean;
ss = sum(sum(ya.*ya));
df = (size(yy,1)-1)*(size(yy,2)-1);
ms_resid = ss/df;
r = diag(repmat(ms_resid,size(yy,2),1));
