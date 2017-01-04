function [pnew, rnew] = Push(AId, Id, p, r, alpha)

pnew = p;
rnew = r;

pnew(Id) = p(Id) + alpha*r(Id);
rnew(Id) = (1 - alpha) * (r(Id) / 2);

GId = find(AId > 0);
Du = sum(AId);
for i = 1 : length(GId)
    tid = GId(i);
    if(tid == Id)
        continue;
    end
    rnew(tid) = r(tid) + ((1 - alpha)/2) * r(Id) * (AId(tid) / Du);
end