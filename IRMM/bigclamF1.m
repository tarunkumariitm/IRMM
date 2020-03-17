function [F] = bigclamF1(clusterLabels, groundTruth)
    n = numel(clusterLabels);
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    C_star_count = max(groundTruth);
    C_cap_count = max(clusterLabels);
    
    if (min(groundTruth ==0))
    groundTruth = groundTruth + 1;
    end
    
    if (min(clusterLabels ==0))
    clusterLabels = clusterLabels + 1 ;
    end
    
    F1_gt_array = zeros(C_star_count, 1);
    F1_cluster_array = zeros(C_cap_count, 1);
    
    for c_star = 1 : C_star_count
        for c_cap = 1: C_cap_count
            TP = 0;
            TN = 0;
            FP = 0;
            FN = 0;
            F1 = 0;
            for i = 1 : n
                for j = i + 1 : n
                    if ((groundTruth(i) == groundTruth(j)) && groundTruth(i) == c_star)
                        if ((clusterLabels(i) == clusterLabels(j)) && clusterLabels(i) == c_cap)
                            TP = TP + 1;
                        else
                            FN = FN + 1;
                        end
                    elseif ((groundTruth(i) ~= groundTruth(j)) && groundTruth(i) == c_star)
                        if ((clusterLabels(i) == clusterLabels(j)) && clusterLabels(i) == c_cap)
                            FP = FP + 1;
                        else
                            TN = TN + 1;
                        end
                    end
                end
            end
            Precision_1 = TP/(TP + FP);
            Recall_1 = TP/(TP + FN);
            F1 = 2 * (Precision_1 * Recall_1 ) / (Precision_1 + Recall_1);
            if (F1 > F1_gt_array(c_star))
                F1_gt_array(c_star) = F1;
            end
        end
    end
    
    
    
    
    for c_cap = 1 : C_cap_count
        for c_star = 1: C_star_count        
            TP = 0;
            TN = 0;
            FP = 0;
            FN = 0;
            F1 = 0;
            for i = 1 : n
                for j = i + 1 : n
                    if ((clusterLabels(i) == clusterLabels(j)) && clusterLabels(i) == c_cap)
                        if ((groundTruth(i) == groundTruth(j)) && groundTruth(i) == c_star)
                            TP = TP + 1;
                        else
                            FN = FN + 1;
                        end
                    elseif ((clusterLabels(i) ~= clusterLabels(j)) && clusterLabels(i) == c_cap)
                        if ((groundTruth(i) == groundTruth(j)) && groundTruth(i) == c_star)
                            FP = FP + 1;
                        else
                            TN = TN + 1;
                        end
                    end
                end
            end
            Precision_1 = TP/(TP + FP);
            Recall_1 = TP/(TP + FN);
            F1 = 2 * (Precision_1 * Recall_1 ) / (Precision_1 + Recall_1);
            if (F1 > F1_cluster_array(c_cap))
                F1_cluster_array(c_cap) = F1;
            end 
        end
    end
    F = (0.5 *(sum(F1_gt_array)/C_star_count)) + (0.5 *(sum(F1_cluster_array)/ C_cap_count));
    
end
