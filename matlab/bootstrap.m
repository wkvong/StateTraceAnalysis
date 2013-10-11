function yb = bootstrap (y)
% draws bootstrap sample from data
if iscell(y) % y = cell array of NSUBJ x NCOND matrices
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
