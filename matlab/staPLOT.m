function staPLOT (data, model, groups, labels, axislabels, axislimits)
% generates a state-trace plot


if ~iscell(data)
    %{
    if sum(sum(isnan(data)))>0 % if data contain NaNs then it's in Henson format
        data = hen2gen (data); % convert Henson format to "general" format
        ys = outSTATS(data);
    end
    %}
    ys = outSTATS (data); % assumes general format
elseif isstruct(data{1})
    ys = data; % if structured then already in stats form
else
    ys = staSTATS(data); % otherwise assume within-subjects data and get stats
end

x = ys{1}.means; 
y = ys{2}.means; 

if ismatrix(ys{1}.n)
%    cx = sqrt(diag(ys{1}.lm)./diag(ys{1}.n));
%    cy = sqrt(diag(ys{2}.lm)./diag(ys{2}.n));
    cx = sqrt(diag(ys{1}.cov)./diag(ys{1}.n)); % between-subjects error bars
    cy = sqrt(diag(ys{2}.cov)./diag(ys{2}.n));
elseif isvector(ys{1}.n)
    cx = sqrt(diag(ys{1}.lm)./ys{1}.n); % within-subjects error bars
    cy = sqrt(diag(ys{2}.lm)./ys{2}.n);
else
    cx = sqrt(diag(ys{1}.lm)/ys{1}.n);
    cy = sqrt(diag(ys{2}.lm)/ys{2}.n);
end
if nargin < 3 || isempty(groups)
    groups = {1:numel(ys{1}.means)};
end
if nargin < 2
    model = [];
end
if nargin < 4 || isempty(labels)
    labels=cell(1,numel(groups));
    for i=1:numel(labels)
        labels{i} = ['Condition ' num2str(i)];
    end
end
if nargin < 5 || isempty(axislabels)
    axislabels = cell(1,2);
    for i=1:numel(axislabels)
        axislabels{i} = ['Outcome Variable ' num2str(i)];
    end
end
if nargin < 6
    axislimits = {};
end
   
% plot data
plotdata (x, y, groups, 0);

% plot error bars
errorbarx (x, y, cx, cy);
errorbary (x, y, cx, cy);

% plot model
if ~isempty(model)
    xm = model{1}; ym = model{2};
    ix = tiesort(xm,ym);
    plot (xm(ix), ym(ix), 'k:'); 
end

% re-plot data over error bars
plotdata (x, y, groups, 1);
hold off

% set axis limits
if ~isempty(axislimits)
    xlim(axislimits{1}); ylim(axislimits{2});
else
    xlim('auto'); ylim('auto');
end

% output axis labels
a = 'xlabel('; a = [a '''' axislabels{1} '''' ');']; eval(a);
a = 'ylabel('; a = [a '''' axislabels{2} '''' ');']; eval(a);

% output legend
a = 'legend(';
for i=1:numel(labels)
    a=[a '''' labels{i} '''' ','];
end
a=[a '''location'', ''southeast'');'];
eval(a);

function plotdata(x,y,groups,flag)
msize=10;
for igroup = 1:numel(groups)
    a = x(groups{igroup});
    b = y(groups{igroup});
    if igroup==1
        h=plot (a, b, 'ko', 'markerfacecolor', 'k', 'markersize', msize); hold on
    elseif igroup==2
        h=plot (a, b, 'ko', 'markerfacecolor', 'w', 'markersize', msize);
    elseif igroup==3
        h=plot (a, b, 'k^', 'markerfacecolor', 'k', 'markersize', msize);
    elseif igroup==4
        h=plot (a, b, 'k^', 'markerfacecolor', 'w', 'markersize', msize);
    elseif igroup==5
        h=plot (a, b, 'ks', 'markerfacecolor', 'k', 'markersize', msize);
    else 
        h=plot (a, b, 'ks', 'markerfacecolor', 'w', 'markersize', msize);
    end
    if flag
        DeleteLegendEntry (h);
    end 
end

function errorbarx (x, y, cx, cy)
% plots errorbars in x direction
yrange=max(y)-min(y)+3*max(cy); 
ticklength=yrange/20;
for i=1:numel(x)
    a(1) = x(i)-cx(i);
    a(2) = x(i)+cx(i);
    b(1) = y(i); 
    b(2) = y(i);
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
% plot ticks
    for j=1:2
        aa(1)=a(j); 
        aa(2)=a(j);
        bb(1)= b(j)-ticklength/2;
        bb(2)= b(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end

function errorbary (x, y, cx, cy)
% plots error bars in y direction
xrange=max(x)-min(x)+3*max(cx); 
ticklength=xrange/20;
for i=1:numel(x)
    a(1) = x(i);
    a(2) = x(i);
    b(1) = y(i)-cy(i); 
    b(2) = y(i)+cy(i);
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
% plot ticks
    for j=1:2
        bb(1)=b(j); 
        bb(2)=b(j);
        aa(1)= a(j)-ticklength/2;
        aa(2)= a(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end

function DeleteLegendEntry (h)
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%{
function newdata = hen2gen (subjectdata, subno)
% converts Henson data to general format
% subno is optional subject number to be inserted
if nargin==1
    subno=1;
end
newdata=[];
for ivar = 1:2
    for icond = 1:3
        j = (ivar-1)*3+icond;
        k = find(~isnan(subjectdata(:,j)));
        u = subjectdata(k,j);
        n = ones(size(u,1),1);
        a = [n*subno n*icond n*ivar u];
        newdata = [newdata; a];
    end
end
%}

function ix = tiesort (xx, yy)
% sorts y values in increasing x-order
% within blocks of tied x-values, sorts y values in increasing y-order
bignumber=100;
x=round(xx*bignumber)/bignumber;
y=round(yy*bignumber)/bignumber;
t=1:numel(x);
z=[x y t'];
a=sortrows(z,[1 2]);
ix=a(:,3);

