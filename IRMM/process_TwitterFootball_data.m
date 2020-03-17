groundTruth = load('./TwitterFootball/classLabels');
removeRows = [22    67    70   104   113   115   117   118   147   163   172   178   192   239];
groundTruth = removerows(groundTruth, removeRows);
n = length(groundTruth);
c = max(groundTruth);


listHypergraph = load('./TwitterFootball/listsMergedHypergraph.mtx');
tempList = spconvert(listHypergraph);
temp = zeros(size(tempList, 1), size(tempList, 2));
temp(find(tempList ~= 0)) = 1;
H{1} = temp;


H{1} = removerows(H{1}, removeRows);
%%%H{2} = removerows(H{2}, removeRows);

remove_H1 = [    2154        3060        3062        3069        3101        3123        3133        3369  3442        3444        3447        3449        3450        3451];
temp = H{1}';
temp = removerows(temp, remove_H1);
H{1} = temp';

contentData = [H{1}];
G = [];
%%%contentData = [H{1} H{2}];

