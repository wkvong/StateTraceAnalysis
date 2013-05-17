function ys = outSTATS (data)
% return statistics for data in "general format"
%subject = data(:,1);
cond = data(:,2); ucond = unique(cond);
var = data(:,3); uvar=unique(var);
within = data(:,4:end); 
for ivar = 1:numel(uvar)
    ys{ivar}.means=[]; 
    ys{ivar}.cov=[]; 
    ys{ivar}.covs=[]; 
    ys{ivar}.n=[]; 
    ys{ivar}.weights=[];
    a=[]; b=[]; c=[];
    for icond = 1:numel(ucond)
        k = find(var==uvar(ivar) & cond==ucond(icond));
        x = within(k,:);
        u = staSTATS(x);
        ys{ivar}.means = [ys{ivar}.means; u.means];
        a{icond} = u.cov;
        b{icond} = u.covs;
        c{icond} = repmat(u.n, size(u.cov));
    end
    s='ys{ivar}.cov = blkdiag(';
    for i=1:numel(ucond)
        s=[s 'a{' num2str(i) '}']; 
        if i < numel(ucond)
            s=[s ','];
        else
            s=[s ');'];
        end
    end
    eval(s);
    s='ys{ivar}.covs = blkdiag(';
    for i=1:numel(ucond)
        s=[s 'b{' num2str(i) '}'];  
        if i < numel(ucond)
            s=[s ','];
        else
            s=[s ');'];
        end
    end
    eval(s);   
    s='ys{ivar}.n = blkdiag(';
    for i=1:numel(ucond)
        s=[s 'c{' num2str(i) '}']; 
        if i < numel(ucond)
            s=[s ','];
        else
            s=[s ');'];
        end
    end
    eval(s); 
    ys{ivar}.weights = diag(ys{ivar}.n)./diag(ys{ivar}.cov);
end
