function [x fit g2fit] = staMRMN (y, E, flag)
% fits partial order model to multinomial data
% y = matrix of counts (rows=conditions, columns=response categories)
% E is the partial order specified as either (a) a cell array, (b)
% adjacency matrix, or (c) inequality matrix
% if flag==1 then do ML optimisation

if nargin==2
    flag=0;
end
nsum = repmat(sum(y,2),1,size(y,2));
p = y./nsum; 
p=reshape(p,numel(p),1);
weights = diag(reshape(nsum,numel(nsum),1));
%weights=eye(numel(nsum));

if iscell(E)
    % construct adjacency matrix
    adj = cell2adj(1:size(y,1), E);
    % construct inequality matrix for first column
    a = adj2ineq(adj);
elseif ~ismember(-1,E) % check if not inequality matrix
    a = adj2ineq(E);
else
    a = E; %E is an inequality matrix already
end

% augment for all columns
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
%    options = optimset ('Algorithm', 'interior-point', 'display','off', 'TolFun', 1e-4, 'TolX', 1e-6, 'UseParallel','always');
    warningstate = warning ('off', 'all'); % turn off rank deficient matrix warnings
    x = fmincon(myfun, x1, A, b, Aeq, beq, lb, ub, [], options);
    warning (warningstate); %reinstate warnings
    u = C*x-d; fit=u'*u;
end
x = reshape(x,size(y,1),size(y,2)).*nsum; % convert to matrix of counts
g2fit = G2(x, y);

function g = G2(x, y)
v = x;
k = find(x == 0); [i j]=ind2sub(size(x),k); yc = sum(y,2)/sum(sum(y));
for l=1:numel(k)
    v(k(l)) = yc(i(l));
end
k=find(y > 0); u = y; u(k) = u(k).*log(y(k)./v(k));
g = sum(sum(u))*2;

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


