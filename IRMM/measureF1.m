function [F] = measureF1(clusterLabels, groundTruth)
    n = numel(clusterLabels);
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    C_star = max(groundTruth);
    C_cap = max(clusterLabels);
    
    for i = 1 : n
        for j = i + 1 : n
            if (clusterLabels(i) == clusterLabels(j))
                if (groundTruth(i) == groundTruth(j))
                    TP = TP + 1;
                else
                    FP = FP + 1;
                end
            else
                if (groundTruth(i) == groundTruth(j))
                    FN = FN + 1;
                else
                    TN = TN + 1;
                end
            end
        end
    end
    
    Precision_1 = TP/(TP + FP);
    Recall_1 = TP/(TP + FN);
    F_1 = 2 * (Precision_1 * Recall_1 ) / (Precision_1 + Recall_1);
    
    
    
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    for i = 1 : n
        for j = i + 1 : n
            if (groundTruth(i) == groundTruth(j))
                if (clusterLabels(i) == clusterLabels(j))
                    TP = TP + 1;
                else
                    FP = FP + 1;
                end
            else
                if (clusterLabels(i) == clusterLabels(j))
                    FN = FN + 1;
                else
                    TN = TN + 1;
                end
            end
        end
    end

    Precision_2 = TP/(TP + FP);
    Recall_2 = TP/(TP + FN);
    F_2 = 2 * (Precision_2 * Recall_2 ) / (Precision_2 + Recall_2);
    F = (F_1 + F_2) /2;
    
    
    %RI = 2 * (TP + TN) / (n * (n - 1));
end
