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
