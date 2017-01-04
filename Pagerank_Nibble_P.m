function [Cluster, Conductance] = Pagerank_Nibble_P(A, NodeId, conductance_threshold)

% vol0 \in [Vol/2, Vol]

N = length(A);
Vol = sum(sum(A));
alpha = (conductance_threshold^2) / (225*log(100*sqrt(Vol)));%(conductance_threshold^2) / (log(Vol));
B = ceil(log(Vol)/log(2));
epsilon = (1/(48*B))*(2^(-2));%1/(10*log(vol0));1/(10*Vol);%
c = 1/8;

Deg = sum(A)';

p = ApproximatePagerank(A, NodeId, alpha, epsilon);

pt = p ./ Deg;
%ptt = pt.*(pt >= c / vol0);
[ps, Id] = sort(pt, 'descend');
Num_gep0 = sum(ps>0);
%Num_gep0 = min(20, Num_gep0);

Degs = Deg(Id);
Conductance = [];
for i = 1 : Num_gep0
    Conductance = [Conductance (sum(Degs(1:i)) - sum(sum(A(Id(1:i), Id(1:i))))) / min(sum(Degs(1:i)), Vol - sum(Degs(1:i)))];
end
% largest cluster
% [MaxCond, MaxCondId] = find(Conductance >= 1 - conductance_threshold);
% if(length(MaxCondId) == 0)
%     Cluster = [];
% else
%     Cluster = Id(1:MaxCondId(end));
% end

% best conductance
[MinCond, MinCondId] = min(Conductance);
Cluster = Id(1:MinCondId);
Conductance = MinCond;


