function path = FindingPath(mov, DistMat, Index, Threshold)
%**************************************************************************
% Construct the Graph using Threshold and MinSpanningTree
[nRows nColumns nColors nFrames] = size(mov);
tempDistMat = DistMat;
%Threshold = mean(DistMat(:));
%Threshold = 0.3;
tempDistMat(tempDistMat>Threshold) = 0;

MST = minspantree(graph(DistMat));
[startInd, endInd] = findedge(MST);

mstDistMat = zeros(nFrames,nFrames);
for i=1:nFrames-1
    mstDistMat(startInd(i),endInd(i)) = DistMat(startInd(i),endInd(i));
end

for i=1:nFrames
    for j=1:nFrames
        if (mstDistMat(i,j)~=0)
            tempDistMat(i,j) = mstDistMat(i,j);
            tempDistMat(j,i) = mstDistMat(i,j);
        end
    end
end

NewDistMat = tempDistMat;
NewDistMat = sparse(NewDistMat);
% newGraph = biograph(NewDistMat);
% view(newGraph);
[dist, path, pred] = graphshortestpath(NewDistMat,Index);
% Save the Path
end
%**************************************************************************

