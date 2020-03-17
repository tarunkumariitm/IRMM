clear all;
maxIter = 1;
num_K = 1;

current_dataset = 'Football';

method = ''; %set to 'clique' for Clique-Louvain method 


process_TwitterFootball_data; % This is a multiview data. This call will create all the layers/views
numHGraphs = length(H); % Count of Hypergraph layers
numGGraphs = length(G); % Count of graph layers
numGraphs = numHGraphs + numGGraphs; % Total number of layers
alpha = [1 0 0 0 0]; % Non zero entry indicates the layer in use. We used the first layer (a hypergraph).K = 20;

maxIter = 12; % set 1 for Hypergraph Louvain. 
K_list = [2, 4, 8, 16, 18, 20, 32]; % varying number of clusters
num_K = length(K_list);
if strcmp(method, 'clique')
    maxIter = 1;
end

counter = 0;

cluster_over_k_iterations = zeros(n, num_K);
for k_iter = 1:num_K
K = K_list(k_iter)
counter = counter + 1;
    numHGraphs = length(H);
    numGGraphs = length(G);
    numGraphs = numHGraphs + numGGraphs;
    n = size(H{1}, 1);

    numHyperedges = 0;
    for i = 1 : numHGraphs
        numHyperedges = numHyperedges + size(H{i}, 2);
    end
    W = ones(numHyperedges, 1);

    H_t = [];
    for i = 1 : size(H, 2)
        H_t = [H_t H{i}];
    end

    La = {};
    invDv = {};
    invDe = {};
    seenHyperedges = 0;
    combinedAdjacency = zeros(n);
    start = cputime;
     flag = 1;
    for i = 1:numHGraphs
        [La{i}, invDv{i}, invDe{i}] = HypergraphUtils.computeAdjacencyMatrix_for_modularity_hypergraph(H{i}, W(seenHyperedges + 1 : seenHyperedges + size(H{i}, 2)), 1);
        if strcmp(method, 'clique')
            [La{i}] = HypergraphUtils.computeCliqueMatrix_for_modularity_hypergraph(H{i}, W(seenHyperedges + 1 : seenHyperedges + size(H{i}, 2)), 1);
        end
         seenHyperedges = seenHyperedges + size(H{i}, 2);
    end

    for i = 1:numGGraphs
        La{numHGraphs+i} = G{i};
    end

    F1 = zeros(1, maxIter);
    bigF1 = zeros(1, maxIter);

    W_all = zeros(numHyperedges, maxIter);
    clusterAssignments = zeros(n, maxIter);
    L_g = zeros(n);
    for i = 1 : numGGraphs
        L_g = L_g + alpha(i + numHGraphs) * La{numHGraphs + i};
    end
    L = L_g;

    for i = 1 : numHGraphs
        L = L + alpha(i) * La{i};
    end



    bin_exp = zeros(11,maxIter);
    for iternum = 1 : maxIter
        iternum
        [COMTY, ending] = cluster_jl_cpp(L, 1,1,0,0);
        cluster_vector = COMTY.COM{1}';
        clusters = cluster_vector;
        if (iternum == maxIter)
            max1 = max(cluster_vector)
            clusters = cluster_agglomorative_hypergraph(K, cluster_vector,L);
        end
        numClusters = max(clusters)
        clusterAssignments(:, iternum) = clusters;
        cluster_over_k_iterations(:, counter) = clusters;
        F1(iternum) = measureF1(clusters, groundTruth);
        bigF1(iternum) = bigclamF1(clusters, groundTruth);
        [newW, bins] = HypergraphUtils.updateWeights(H_t, W, clusters);
        sum(bins);
        bin_exp(:,iternum) = bins;
        
        W = newW;
        W_all(:, iternum) = newW;
        seenHyperedges = 0;
        for i = 1:numHGraphs
            [La{i}] = HypergraphUtils.computeIncrementalAdjacency_hypergraph(H{i}, W(seenHyperedges + 1 : seenHyperedges + size(H{i}, 2)), invDv{i}, invDe{i});
            seenHyperedges = seenHyperedges + size(H{i}, 2);
        end

        L = L_g;
        for i = 1 : numHGraphs
            L = L + alpha(i) * La{i};
        end
    end
    stop = cputime;
    cputime_elapsed = stop - start

    weightChange = zeros(1, iternum);
    weightChange(1) = dist(ones(1, numHyperedges), W_all(:, 1));
    for i = 2 : iternum
        weightChange(i) = dist(W_all(:, i - 1)', W_all(:, i));
    end

    assignChange = zeros(1, iternum);
    for i = 2 : iternum
        assignChange(i) = measureRandIndex(W_all(:, i), W_all(:, i - 1));
    end
    minMaxDiff = max(W_all) - min(W_all);


F1
bigF1
end
% 
% csvwrite([folderName, 'RI.csv'], RI);
% 
% csvwrite([folderName, 'purity.csv'], averagePurity);
% csvwrite([folderName, 'W_all.csv'], W_all);
% csvwrite([folderName, 'clusters.csv'], clusterAssignments);
% csvwrite([folderName, 'weightChange.csv'], weightChange);
% csvwrite([folderName, 'assignChange.csv'], assignChange);
% csvwrite([folderName, 'minMaxDiff.csv'], minMaxDiff);
% csvwrite([folderName, 'num_of_points.csv'], n);
% csvwrite([folderName, 'num_of_clusters.csv'], numClusters);
% csvwrite([folderName, 'alpha.csv'], alpha);
% 
% sorted_cluster_bucket = sort(cluster_bucket);
% histogram(sorted_cluster_bucket)
% 
