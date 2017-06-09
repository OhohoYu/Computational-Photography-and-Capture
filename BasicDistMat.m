function [DistMat, ProbMat] = BasicDistMat(mov)
%**************************************************************************
% Compute Basic Distance Matrix and Probability Matrix:
[nRows nColumns nColors nFrames] = size(mov);
DistMat = zeros(nFrames,nFrames);
ProbMat = zeros(nFrames,nFrames);
for j = 1:nFrames
    for i = 1:nFrames
        difFrame = mov(:,:,:,i) - mov(:,:,:,j);
        distFrame = abs(difFrame);
        distValue = sum(distFrame(:))/(nRows*nColumns);
        DistMat(i,j) = distValue;
    end
end
DistMat = DistMat/max(max(DistMat));
sigma = sum(DistMat(:))/sum(sum(DistMat~=0));
ProbMat(1:end-1,:) = exp(-DistMat(2:end,:)/sigma);
sumVec = sum(ProbMat(1:end-1,:),2);
ProbMat(1:end-1,:) = ProbMat(1:end-1,:)./repmat(sumVec,1,nFrames);
end
%**************************************************************************

