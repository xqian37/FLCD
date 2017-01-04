function [ClusterTT] = PageRank_Nibble_Bath(A, Cond_threshold, Overlapping_threshold)

A = sparse(A);


N = length(A);

k = 1;
Conductance = [];
for i = 1 : N
    [ClusterT, CondT] = Pagerank_Nibble_P(A, i, Cond_threshold);
    if(length(ClusterT) <= 2 || CondT > Cond_threshold)
        continue;
    end
    Cluster{k} = ClusterT;
    Conductance = [Conductance CondT];
    %disp(['Cluster: ' num2str(k) ' length: ' num2str(length(ClusterT)) ' i ' num2str(i) ' Cond ' num2str(CondT)]);
    k = k + 1;
end


% remove highly overlapped clusters
NumC = length(Cluster);
for i = 1 : length(Cluster)
    len1 = length(Cluster{i});
    if(len1 <= 2)
        continue;
    end
    for j = i+1 : length(Cluster)
        len2 = length(Cluster{j});
        temp = intersect(Cluster{i},Cluster{j});
        overlap = length(temp);
        overlap_score = overlap^2 / (len1*len2);
        if(overlap_score >= Overlapping_threshold)
            Cluster{j} = [];
            NumC = NumC-1;
            %disp(['Reomve Cluster ' num2str(j) ' duplicate with Cluster ' num2str(i)]);
        end
    end
end

k = 1;
for i = 1 : length(Cluster)
    if(length(Cluster{i})<=2)
        continue;
    end
    ClusterTT{k} = Cluster{i};
    k = k + 1;
end







