function [ClusterInd] = PutBackIn(SubInd, Dense, A)
%
T1 = (2+1/length(SubInd))*(Dense/2);
T2 = Dense/2;

NodeSet = [];
for i = 1 : length(SubInd)
    Temp = find(A(SubInd(i),:)==1);
    NodeSet = [NodeSet Temp];
end

UNSet = unique(NodeSet);
AddIn = [];
for i = 1:length(UNSet)
    Dt = sum(NodeSet==UNSet(i));
    if(Dt >= T1 && ~any(SubInd==UNSet(i)))
        AddIn = [AddIn; UNSet(i)];
    end
end

ClusterInd = [SubInd; AddIn];