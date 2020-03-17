classdef HypergraphUtils
    methods (Static)
        function [Pi_mix, Theta_mix] = combineRandomWalk(pi, theta, alpha)
            % pi is nXd matrix where n is the number of instances and d is
            % the number of views

            if length(alpha)==1
                Pi_mix = pi;
                Theta_mix = theta{1};
                return;
            end

            n = size(pi, 1);
            Pi_mix = sparse(n, 1);

            for i=1:length(alpha)
                Pi_mix = Pi_mix + alpha(i)*pi(:, i);
            end

            beta = zeros(size(pi));

            for i = 1:length(alpha)
                beta(:, i) = alpha(i)*(pi(:, i)./Pi_mix);
            end

            clear pi;

            Theta_mix = sparse(n, n);
            for i=1:length(alpha)
                Theta_mix = Theta_mix + diagonalize(beta(:, i))*theta{i};
            end

        end
        
        %% stochastic matrix for hypergraph
        function [theta, Pi] = computeStochaticMatrix_hypergraph(H, W)    
            
            dv = sum(H * diag(W), 2); %vertex degree vector
            invdv = dv.^(-1);
            Pi = dv/sum(dv);
            clear dv;
            invDv = diagonalize(invdv);
            clear invdv;
            de = sum(H, 1)-1; %edge degree vector
            De = diagonalize(de); %edge degree diagonal matrix
            clear de;

            theta = invDv*H*De^(-1)*H';  % stochastic matrix
            theta(speye(size(theta))==1) = 0;        
        end
        
         %% stochastic matrix for graph
        function [theta, Pi] = computeStochasticMatrix_graph(G)
            %G is an adjacency matrix
            outDegree = sum(G, 2);
            invOutDegree = outDegree.^(-1);
            Dv = diagonalize(outDegree);
            invDv = diagonalize(invOutDegree);
            %beta = 0.9;
            %n = size(G, 1);

            %theta = beta*invDv*G + (1-beta)*ones(n,n)*(1/n);
            %Pi = beta*Dv/trace(Dv) + (1-beta)*speye(n,n)*(1/n);
            theta = invDv*G;
            theta = full(theta);
            Pi = outDegree/trace(Dv);
            %theta = (Pi^0.5*theta*Pi^-0.5 + Pi^-0.5*theta'*Pi^0.5)/2;
        end
        
        function [L, D] = computeUnnormalizedLaplacian_hypergraph(H, W)
            
        end
            
        function [L, L1, Dv, invDv, invDe] = computeNormalizedLaplacian_hypergraph(H, W)
            n = size(H, 1);
%             W = eye(size(H,2));
            dv = sum(H, 2); %vertex degree vector
            invdv = full(dv).^(-1/2);
            invDv = diag(invdv);
            Dv = diag(dv);
%             clear dv;
%             invDv = full(Dv)^(-1/2);
%             invDv = diagonalize(invdv);
            clear invdv;
            de = sum(H, 1)-1;%edge degree vector
            invde = full(de).^(-1);
            invDe = diag(invde);
%             De = diagonalize(de);
%             invDe = full(De)^(-1);
            clear invde;
%             invDe = diagonalize(invde); %edge degree diagonal matrix
            L1 = invDv * H * diag(W) * invDe * H' * invDv;
            L = eye(n) - L1;
        end
        
        function [L] = computeNormalizedLaplacian_clique(H, W)
            n = size(H, 1);
            A = H * diag(W) * H';
            A(logical(eye(size(A)))) = 0;
            dv = sum(A, 2);
            invdv = full(dv).^(-1/2);
            invDv = diag(invdv);
            L = eye(n) - (invDv * A * invDv);
            %L1 = invDv * H * diag(W) * invDe * H' * invDv;
            %L = eye(n) - L1;
        end

        
        
        function [L1] = computeIncrementalLaplacian_hypergraph(H, W, invDv, invDe)
            L1 = invDv * H * diag(W) * invDe * H' * invDv;
        end
        
        function [L1] = computeIncrementalAdjacency_hypergraph(H, W, invDv, invDe)
            L1 = H * diag(W) * invDe * H';
        end

        
        function [L, L1, D] = computeNormalizedLaplacian_graph(G)
            % G is the weighted adjacency matrix
           n = size(G, 1);
           dv = sum(G, 2); %vertex degree vector
           invdv = full(dv).^(-1/2);
           invdv(find(isinf(invdv))) = 0;
           invDv = diag(invdv);
           L1 = invDv * G * invDv;
           L = eye(n) - L1;
           D = diag(dv);
        end
        
        
        function [A] = computeCliqueMatrix_for_modularity_hypergraph(H, W, weighted);
            %n = size(H, 1);
            %De = sum(H, 1);
            %De = De';
            A = H * diag(W) * H';
            A(logical(eye(size(A)))) = 0;

        end
        function [A, invDv, invDe] = computeAdjacencyMatrix_hypergraph(H, W, weighted)
            % If weighted = 1, use weighted clique formulation
            % If not, use normal
           	n = size(H, 1);
            De = sum(H, 1);
            De = De';
%             if (weighted ~= 0)
%                 W = W .* (1 ./ (De - 1));
% %                 W = W1;
%             end
             de = sum(H, 1)-1; %edge degree vector
             invde = full(de).^(-1);
             invDe = diag(invde);
%             
            dv = sum(H, 2); %vertex degree vector
            invdv = full(dv).^(-1/2);
            invDv = diag(invdv);
            
            A = H * diag(W) * invDe * H';
            % Set diagonal elements to 0
            A(logical(eye(size(A)))) = 0;          
            
        end
        
        
        function [A, invDv, invDe, Dv, norm_W] = computeAdjacencyMatrix_for_modularity_hypergraph(H, W, weighted)
		    % If weighted = 1, use weighted clique formulation
		    % If not, use normal
		   	n = size(H, 1);
		    De = sum(H, 1);
		    De = De';
%		    if (weighted ~= 0)
		        W = W .* (1 ./ (De - 1));
	%                 W = W1;
%		    end
		    norm_W = W;
		    
		    dv = sum(H, 2); %vertex degree vector
		    invdv = full(dv).^(-1/2);
		    invDv = diag(invdv);
		    de = sum(H, 1)-1; %edge degree vector
		    invde = full(de).^(-1);
		    invDe = diag(invde);
		    Dv = diag(dv);
		    
		    A = H * diag(W) * H';
                    %A = A - Dv;  %% This is used in Harini's thesis but not coded
                    % Set diagonal elements to 0
		    A(logical(eye(size(A)))) = 0;
        end

        
        
        function [M] = computeModularityMatrix_hypergraph(H, W, weighted)
            A = HypergraphUtils.computeAdjacencyMatrix_hypergraph(H, W, weighted);
            
            M = HypergraphUtils.computeModularityforAdjacency(A);
            
        end
        
        function [M] = computeModularityMatrix_graph(A)
            % Assume A is weighted adjacency matrix
            M = HypergraphUtils.computeModularityforAdjacency(A);
        end
        
        function [M] = computeModularityforAdjacency(A)
            k = sum(A, 2);
            m = sum(k);
            P = (k * k') / m;
            M = A - P;
        end
        
        function [S, groups, tempGroups] = assignGroups(Y)
            % Assume columns of Y are Eigen vectors in order of decreasing
            % Eigen values
            n = size(Y, 1);
            seen = zeros(1, n);
            [~, tempGroups] = max(Y');
            col = 1;
            while (sum(seen) < n && col <= size(Y, 2))
                for i = 1 : n
                    if (Y(i, col) > 0 && seen(i) == 0)
                        S(i, col) = 1;
                        seen(i) = 1;
                    end
                end
                col = col + 1;
            end
            % If there is an unassigned point, assign it to 1
            for i = 1 : n
                if (sum(S(i,:)) == 0)
                    S(i, 1) = 1;
                end
            end
            % If there is a group with no elements, remove it
            S( :, ~any(S,1) ) = [];
            [~, groups] = max(S');
        end
        
        function [H] = convertAdjacencyToIncidence(A)
            % Assuming unweighted adjacency
            n = size(A, 1);
            H = [];
            e = 1;
            for i = 1 : n
                for j = i + 1 : n
                    if (A(i, j) == 1)
                        H(:, e) = zeros(n, 1);
                        H(i, e) = 1;
                        H(j, e) = 1;
                        e = e + 1;
                    end
                end
            end
        end
        
        function [newW, bins] = updateWeights(H, W, clusters)
            % Assume W is just a vector of weights, not a diagonal matrix
            m = size(H, 2);
            n = size(H, 1);
            k = max(clusters);
            changeInWeights = zeros(m, 1);
            bins = zeros(11,1);
            for i = 1 : m
                nodesInHyperedge = find(H(:, i) == 1);
                nodeClusters = clusters(nodesInHyperedge);
                numNodesInHyperedge = numel(nodesInHyperedge);
                fractions = zeros(k,1);
                for j = 1 : k
                    numNodesInCluster = sum(nodeClusters == j);
%                     if (numNodesInCluster ~= 0)
%                     changeInWeights(i) = changeInWeights(i) + numNodesInCluster * (numNodesInCluster - 1);
                        changeInWeights(i) = changeInWeights(i) + (1 / (numNodesInCluster + 1));
                        fractions(j) = numNodesInCluster/numNodesInHyperedge;%Tarun
                end
                
                fraction = max(fractions);

                if fraction == 1
                    bins(11) = bins(11) +1;
                end
                if (fraction > 0) && (fraction<= 0.1)
                    bins(1) = bins(1) +1;
                end
                if fraction > 0.1 && fraction<= 0.2
                    bins(2) = bins(2) +1;
                end
                if fraction > 0.2 && fraction<= 0.3
                    bins(3) = bins(3) +1;
                end
                if fraction > 0.3 && fraction<= 0.4
                    bins(4) = bins(4) +1;
                end
                if fraction > 0.4 && fraction<= 0.5
                    bins(5) = bins(5) +1;
                end
                if fraction > 0.5 && fraction<= 0.6
                    bins(6) = bins(6) +1;
                end
                if fraction > 0.6 && fraction<= 0.7
                    bins(7) = bins(7) +1;
                end
                if fraction > 0.7 && fraction<= 0.8
                    bins(8) = bins(8) +1;
                end
                if fraction > 0.8 && fraction<= 0.9
                    bins(9) = bins(9) +1;
                end
                if fraction > 0.9 && fraction < 1
                    bins(10) = bins(10) +1;
                end
                         
%                     end
%                 changeInWeights(i) = changeInWeights(i) / numNodesInHyperedge * (numNodesInHyperedge - 1);
                changeInWeights(i) = changeInWeights(i) / (1 / (numNodesInHyperedge + k));
            end
            
            changeInWeights = m * changeInWeights / sum(changeInWeights);
            
            %newW = (changeInWeights + W) / 2;
            del_w_wt = 0.5;
            w_wt = 1 - del_w_wt;
            newW = ((del_w_wt .* changeInWeights) + (w_wt .* W));
        end
                
                
        function [H_combined] = combineToIncidenceMatrix(H, G)
            H_combined = [];
            for i = 1 : size(H, 2)
                H_combined = [H_combined H{i}];
            end
            for i = 1 : size(G, 2)
                H_combined = [H_combined HypergraphUtils.convertAdjacencyToIncidence(G{i})];
            end
        end
                    
        function [L_combined] = combineLaplacians(H, G, W, alpha)
            numHGraphs = size(H, 2);
            numGGraphs = size(G, 2);
            n = size(H{1}, 1);
            L_combined = zeros(n);
            
            seenEdges = 0;
            for i = 1 : numHGraphs
                [~, L1, ~] = HypergraphUtils.computeNormalizedLaplacian_hypergraph(H{i}, W(seenEdges + 1 : seenEdges + size(H{i}, 2)));
                seenEdges = seenEdges + size(H{i}, 2);
                L_combined = L_combined + alpha(i) * L1;
            end
            
            
            for i = 1 : numGGraphs
                [~, L1, ~] = HypergraphUtils.computeNormalizedLaplacian_graph(G{i});
                L_combined = L_combined + alpha(i + numHGraphs) * L1;
            end
            
            L_combined = eye(n) - L_combined;
        end
            
            
%         function [] = computeModularityMatrix(H, G, alpha)
%             % H: cell array of hypergraph incidence matrices
%             % G: cell array of graph adjacency matrices (weightes)
%             % alpha: importance of views
%             
%             numHGraphs = size(H, 2);
%             numGHraphs = size(G, 2);
%             numGraphs = numHGraphs + numGGraphs;
%            
%             
%         end

            
            
        
    end
end