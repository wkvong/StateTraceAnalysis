function ys = outSTATS (data)
% return statistics for data in "general format"
% updated 19/6/2013 in light of modifications to staSTATS
% general format is defined as:
% column 1 = subject number (not used directly)
% column 2 = between-subjects condition
% column 3 = dependent variable (between-subjects)
% columns 4 to end = values for each within-subjects condition
% output (ys) is cell array containing statistics for each level of the
% dependent variable
% each componenent of ys is a structured array consisting of:
% means = ncond(between) x ncond(within) matrix of means
% cov = covariance matrix of size ncond(between) x ncond(within)
% n = matrix of size ncond(between) x ncond(within) containing number of
% subjects
% weights = weight vector of length ncond(between) x ncond(within)
%
% subject = data(:,1);
cond = data(:,2); ucond = unique(cond);
var = data(:,3); uvar=unique(var);
within = data(:,4:end); 
for ivar = 1:numel(uvar)
    ys{ivar}.means=[]; 
    ys{ivar}.cov=[]; 
    ys{ivar}.lm=[];
    ys{ivar}.n=[]; 
    ys{ivar}.weights=[];
    a={}; b={}; c={};
    for icond = 1:numel(ucond)
        k = find(var==uvar(ivar) & cond==ucond(icond));
        x = within(k,:);
        u = staSTATS(x);
        ys{ivar}.means = [ys{ivar}.means u.means];
        a{icond} = u.cov;
        b{icond} = u.lm;
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

    s='ys{ivar}.lm = blkdiag(';
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
    a = diag(ys{ivar}.n)./diag(ys{ivar}.cov);
    ys{ivar}.weights = a';
end
