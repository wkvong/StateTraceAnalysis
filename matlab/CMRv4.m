function [xStar Ffits EStar exitflag] = CMRv4(y, E)
% Version 4 of CMR: processes adj matrix directly (a bit faster)
% implements Oleg Burdakov's branch and bound algorithm for solving the
% Coupled Monotonic Regression problem (version 2)
% y is a cell array of structured output from staSTATS (ie means, weights etc) 
% E is a cell array of the starting partial order model: E for Edges
% for example, if we want x3 <= x2 <= x1 and x4 <= x3 then
% E = {[3 2 1] [4 3]}, and so on.
% Fbar is the fit of the starting model -- 'Inf' is a good start
% returns:
% xystar = best fitting values
% Fbar = weighted least squares fit
% Estar = final partial order model 
% exitflag = vector of exit flags for fits for each variable
L = [];
if iscell(E)
    L(1).E = cell2adj (1:numel(y{1}.means), E);
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
        [xPrime fits exitflag] = staMR(y, Eprime); % this is where lsqlin comes in
        Ffit = sum(fits);
        if Ffit < Fbar
            [flag idx] = Feasible(xPrime); %if flag == true, it's feasible, otherwise give me an index
            if flag
                Fbar = Ffit;
                Ffits = fits;
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

function [flag idx] = Feasible(xx)
% returns true if feasible, and the largest inversion if not
% xx is cell array of possible solutions
% this is the biggest difference from CMRv2
flag = 1; idx = []; zero = 1e-010;
x = zeros(numel(xx{1}),numel(xx));
for ivar = 1:numel(xx), x(:,ivar) = xx{ivar}; end

u = nchoosek(1:size(x,1),2);
d = x(u(:,1),:) - x(u(:,2),:);
k = find(abs(d) <= zero); d(k) = 0;  % set to zero any tiny differences
d = sign(d); % convert to signs
s = sum(abs(d),2)-abs(sum(d,2)); % s > 0 if signs of differences not equal 
k = find(s > 0);
if ~isempty(k)
    flag=0;
    idx = u(k(1),:);   % return the indices of first discordant pair
end