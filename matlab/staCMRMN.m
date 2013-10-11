function [xStar Fbar g2fit EStar] = staCMRMN (y, E, mlflag)
% applies CMR algorithm to multinomial data
% y is ncond x nresp matrix of counts
% E is a cell array of the starting partial order model: E for Edges
% for example, if we want x3 <= x2 <= x1 and x4 <= x3 then
% E = {[3 2 1] [4 3]}, and so on.
% conditions are ordered columnwise with respect to y
% mlflag = 1 then do maximum likelihood fit (may not work)
% returns:
% xystar = best fitting values
% Fbar = weighted least squares fit
% g2fit = G-squared fit if calculated (i.e. if mlflag=1)
% EStar = adjacency matrix of latent partial order

if nargin <= 2
    mlflag=0;
end
if nargin == 1
    E={};
end
    
L = [];
if iscell(E)
    L(1).E = cell2adj (1:size(y,1), E);
else
    L(1).E = E;
end
L(1).F = -Inf;
Fbar = Inf; 
EBar = L(1).E;

while numel(L)>0
    Eprime = L(1).E;
    Ffloor = L(1).F;
    L(1) = []; %remove this node
    if Ffloor < Fbar
        [xPrime fit] = staMRMN1(y, Eprime, mlflag); % this is where lsqlin comes in
        Ffit = fit;
        if Ffit < Fbar
            [flag idx] = Feasible(xPrime); %if flag == true, it's feasible, otherwise give me an index
            if flag
                Fbar = Ffit;
                xBar = xPrime;
                EBar = Eprime;
            else
                %create two branches
                L(end+1).E = Eprime; L(end).E(idx(1),idx(2)) = 1; % (i,j) branch
                L(end).F = Ffit;
                L(end+1).E = Eprime; L(end).E(idx(2),idx(1)) = 1; % (j,i) branch
                L(end).F = Ffit;
            end
        end
    end
end
xStar = xBar;
EStar = EBar;
k = find(xStar<0); xStar(k)=0; %round negative values to zero
%g2fit = G2(xStar, y);
g2fit=NaN;

function [flag idx] = Feasible(x)
% returns true if feasible, and the largest inversion if not
% x is a matrix of expected counts
flag = 1; idx = []; zero = 1e-10;
p = cumsum(x,2)./repmat(sum(x,2),1,size(x,2));
u = nchoosek(1:size(p,1),2);
d = p(u(:,1),:) - p(u(:,2),:);
k = find(abs(d) <= zero); d(k) = 0;  % set to zero any tiny differences
d = sign(d); % convert to signs
s = sum(abs(d),2)-abs(sum(d,2)); % s > 0 if signs of differences not equal 
k = find(s > 0);
if ~isempty(k)
    flag=0;
    idx = u(k(1),:);   % return the indices of first discordant pair
end

function [x fit] = staMRMN1 (y, E, flag)
% fits partial order model to multinomial data
% y = matrix of counts (rows=conditions, columns=response categories)
% E is the partial order specified as an adjacency matrix
%
% convert to row proportions
nsum = repmat(sum(y,2),1,size(y,2));
p = y./nsum; 
p=reshape(p,numel(p),1);
weights = diag(reshape(nsum,numel(nsum),1));
% augment for all columns
a = adj2ineq(E);
A = ineqrep(a, size(y,2));
b = zeros(size(A,1),1);
% make sure sum of rows = 1
Aeq = [];
for j=1:size(y,2)
    Aeq = [Aeq eye(size(y,1))];
end
beq = ones(size(Aeq,1),1);
% set up matrices and bounds
C = sqrt(weights);
d = C*p;
x0 = p;
lb=zeros(size(x0));
ub=ones(size(x0));

% set up lsqlin
options = optimset ('LargeScale','off', 'display','off'); % turn off display and use medium algorithm to avoid warning
[x fit] = lsqlin(C, d, A, b, Aeq, beq, lb, ub, x0, options);
x1 = x;

% finish off with ML optimisation
if flag > 0
    myfun=@(x) ML(x, y);
    options = optimset ('Algorithm', 'interior-point', 'display','off', 'TolFun', 1e-4, 'TolX', 1e-6);
    warningstate = warning ('off', 'all'); % turn off rank deficient matrix warnings
    x = fmincon(myfun, x1, A, b, Aeq, beq, lb, ub, [], options);
    warning (warningstate); %reinstate warnings
    u = C*x-d; fit=u'*u;
end
x = reshape(x,size(y,1),size(y,2)).*nsum; % convert to matrix of counts

function m = ML (x, y)
if isvector (x) %called by fmincon assumes x is vector of probabilities
    yy = reshape(y,numel(y),1);
    k = find(x == 0); p = x; p(k)=1/sum(yy);
    u = yy.*log(p);
    m = -sum(u);
else % assumes x is matrix of expected counts
    p = x./repmat(sum(y,2),1,size(y,2));
    k = find(x == 0); p(k)=1/sum(sum(y));
    u = y.*log(p);
    m = -sum(sum(u));
end

function g = G2(x, y)
v = x;
k = find(x == 0); [i j]=ind2sub(size(x),k); 
yc = sum(y,2)/sum(sum(y));
for l=1:numel(k)
    v(k(l)) = yc(i(l));
end
k=find(y > 0); u = y; u(k) = u(k).*log(y(k)./v(k));
g = sum(sum(u))*2;

function adj = cell2adj  (nodes, E)
% converts a partial order model in cell array form to an adjacency matrix suitable for MR
if nargin==1
    E={};
end
if ~iscell(E)
    E={E};
end
n=numel(nodes);
adj=zeros(n,n);
if ~isempty (E)
    for i=1:numel(E)
        if ~isempty(E{i})
            u = nchoosek(E{i},2);
            for j=1:size(u,1)
                k1=find(nodes==u(j,1));
                k2=find(nodes==u(j,2));
                adj(k1,k2)=1;
            end
        end
    end
end

function Aineq = adj2ineq (adj)
% converts an adjacency matrix to a inequality coefficient matrix (A) for
% MR
[i j] = find(sparse(adj));
Aineq = zeros(numel(i), size(adj,1));
for k=1:numel(i)
    Aineq(k,i(k))=1;
    Aineq(k,j(k))=-1;
end

function b = ineqrep(a, n)
b=[];
for i=1:n
    u=[];
    for j=1:n
        if j<=i
            u=[u a];
        else
            u=[u zeros(size(a))];
        end
    end
    b = [b; u];
end


