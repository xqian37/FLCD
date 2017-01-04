function app = ApproximatePagerank(A, NodeId, alpha, epsilon)

N = length(A);

p = zeros(N,1);
r = zeros(N,1);
r(NodeId) = 1;

Deg = sum(A)';

du = sum(A(NodeId,:));
[Maxr, MaxId] = max(r/du);
while(Maxr >= epsilon)
    [pnew, rnew] = Push(A(MaxId,:), MaxId, p, r, alpha);
    p = pnew;
    r = rnew;
    rt = r ./ Deg;
    [Maxr, MaxId] = max(rt);
%     [Rank RankId] = sort(rt, 'descend');
%     Maxr = Rank(2);
%     MaxId = RankId(2);
end

app = p;