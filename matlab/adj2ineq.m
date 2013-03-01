function Aineq = adj2ineq (adj)
% converts an adjacency matrix to a inequality coefficient matrix (A) for
% MR
[i j s] = find(sparse(adj));
Aineq = zeros(numel(i), size(adj,1));
for k=1:numel(i)
    Aineq(k,i(k))=1;
    Aineq(k,j(k))=-1;
end
