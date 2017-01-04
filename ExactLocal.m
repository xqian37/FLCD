function [Clusters] = ExactLocal(A, local_nodes)

Deg = sum(A);
N = length(A);
Vol = sum(sum(A));

conductance_threshold = 0.5;

% set alpha and epsilon for lazy random walk
alpha = (conductance_threshold^2) / (225*log(100*sqrt(Vol)));%(conductance_threshold^2) / (log(Vol));
B = ceil(log(Vol)/log(2));
epsilon = (1/(48*B))*(2^(-3));

NodeList = 1:N;
Node_Deg = [NodeList' Deg'];

changerate = 1e-5;
[MaxDeg MaxDegId] = max(Node_Deg(:,2));
clusterid = 1;
while(length(MaxDeg)>0)
    
    if(MaxDeg < 3)
        Node_Deg(MaxDegId,:) = [];
        [MaxDeg MaxDegId] = max(Node_Deg(:,2));
        continue;
    end
    
    NodeId = Node_Deg(MaxDegId,1);
    p = ApproximatePagerank(A, NodeId, alpha, epsilon);
    ind = find(p>0);
    pt = p(ind);
    [ord indd] = sort(pt,'descend');
    
    if(length(indd)>local_nodes)
        count = local_nodes;
    else
        count = length(indd);
    end
    wind = ind(indd(1:count));
    SubA = A(wind,wind);
    SubD = Deg(wind);
    
    %step 1
    [ClustersT{clusterid} Conductance(clusterid)] = Low_Conductance(SubA, SubD, wind, 0);
    SSA = A(ClustersT{clusterid},ClustersT{clusterid});
    DDT = ones(1,length(ClustersT{clusterid}));
    
    %step 2
    [ClustersTT{clusterid}, Density(clusterid)] = Low_Conductance(SSA, DDT, ClustersT{clusterid}, 0);
    [Clusters{clusterid}] = PutBackIn(ClustersTT{clusterid}, Density(clusterid), A);
    
    clusterid = clusterid + 1;
    
    %post-processing
    if(Density(clusterid-1) < 1000)
        removeid = MaxDegId;
    else
        removeid=MaxDegId;
        for i = 1 : length(Clusters{clusterid-1})
            if(ismember(Clusters{clusterid-1}(i),Node_Deg(:,1))) %&& Conductance(clusterid-1) >= 0.65
                removeid = [removeid find(Node_Deg(:,1)==Clusters{clusterid-1}(i))];
            end
        end
    end
    Node_Deg(removeid,:) = [];
    
    [MaxDeg MaxDegId] = max(Node_Deg(:,2));
    
    disp([num2str(length(Node_Deg)) ' left!'])
    
    
end













