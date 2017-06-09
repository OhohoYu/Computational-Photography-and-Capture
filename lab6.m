close all;
clear all;

% Load Frames
movOriginal = load_sequence_color('./gjbLookAtTargets','gjbLookAtTarget_', 0, 71, 4, 'jpg');
[nRows, nColumns, nColors, nFrames] = size(movOriginal);

for i = 1:nFrames
    mov(:,:,:,i) = imresize(movOriginal(:,:,:,i),0.3);
end
[nRowsDown, nColumnsDown, nColors, nFrames] = size(mov);

% movOriginal = load_sequence_color('./myself','myself_', 0, 19, 3, 'jpg');
% [nRows, nColumns, nColors, nFrames] = size(movOriginal);
% 
% for i = 1:nFrames
%     mov(:,:,:,i) = imresize(movOriginal(:,:,:,i),0.5);
% end
% [nRowsDown, nColumnsDown, nColors, nFrames] = size(mov);


% DistMat = BasicDistMat(mov);
% %DistMat = DP_DistMat(mov);
% Threshold = mean(DistMat(:));
% path = FindingPath(mov, DistMat, 1, Threshold);
% %flow = CreateFlow(mov);
% flow = load('newFlow.mat');
% flow = flow.newFlow;
%**************************************************************************
% Load the Pre-Computed Data:
DistMat = load('DistDown.mat');
DistMat = DistMat.DistMat;
ProbMat = load('ProbDown.mat');
ProbMat = ProbMat.ProbMat;
flow = load('flows.mat');
flow = flow.flows_a;
Threshold = mean(DistMat(:));

% path = load('meanPath.mat');
% path = path.path;

% DistMat = load('DPDistMat.mat');
% DistMat = DistMat.DistMat1;
% ProbMat = load('DPProbMat.mat');
% ProbMat = ProbMat.ProbMat1;
% flow = load('flows.mat');
% flow = flow.flows_a;
% path = load('DPPath.mat');
% path = path.path;
%**************************************************************************
% Single Click
% Run a simple Test Program:
% The 1st click is to define your starting pixel
% while the 2nd click is to locate your destination,
% then the relevant frames would be played in the window.
% while(1)
% imshow(mov(:,:,:,1));
% [tempX, tempY] = ginput(2);
% startPixel = [tempY(1), tempX(1)];
% expectedPixel = [tempY(2), tempX(2)];
% EndPixels = zeros(length(path),2);
% for i=1:length(path)
%     currPath = path{i};
%     endPixel = startPixel;
%     if length(currPath)>1
%         for j=1:length(currPath)-1
%             if currPath(j)>currPath(j+1)
%                 k = (currPath(j)-1)*(currPath(j)-2)/2 + currPath(j+1);
%                 endInd = round(endPixel);
%                 currFlowX = flow(endInd(1),endInd(2),1,k);
%                 currFlowY = flow(endInd(1),endInd(2),2,k);
%                 endPixel = endPixel + [currFlowY, currFlowX];
%             else
%                 k = (currPath(j+1)-1)*(currPath(j+1)-2)/2 + currPath(j);
%                 endInd = round(endPixel);
%                 currFlowX = -flow(endInd(1),endInd(2),1,k);
%                 currFlowY = -flow(endInd(1),endInd(2),2,k);
%                 endPixel = endPixel + [currFlowY, currFlowX];
%             end
%         end
%         EndPixels(i,:) = endPixel;
%     else
%         EndPixels(i,:) = endPixel;
%     end
% end
% 
% DistEndPixels = sqrt(sum((EndPixels-repmat(expectedPixel,length(path),1)).^2,2));
% [minDist, minInd] = min(DistEndPixels);
% FinalPath = path{minInd}
% for i=1:length(FinalPath)
%     imshow(mov(:,:,:,FinalPath(i)));
% end
% 
% % Render Slow Motion Interpolation:
% PathMov = MotionInterp(mov, FinalPath);
% for i=1:length(FinalPath)-1
%     imshow(PathMov(:,:,:,i));
% end
% pause(2);
% end

% *************************************************************************
% Multiple Clicks
while(1)
nClicks = input('Input the Number of Clicks: ');
imshow(mov(:,:,:,1));
[tempX, tempY] = ginput(nClicks); 
path = FindingPath(mov, DistMat, 1, Threshold);
hold on;
plot(tempX, tempY, 'ro-');
hold off;
pause(2);
for count = 1:nClicks-1
EndPixels = zeros(length(path),2);
startPixel = [tempY(count),tempX(count)];
expectedPixel = [tempY(count+1),tempX(count+1)];
for i=1:length(path)
    currPath = path{i};
    endPixel = startPixel;
    if length(currPath)>1
        for j=1:length(currPath)-1
            if currPath(j)>currPath(j+1)
                k = (currPath(j)-1)*(currPath(j)-2)/2 + currPath(j+1);
                endInd = round(endPixel);
                currFlowX = flow(endInd(1),endInd(2),1,k);
                currFlowY = flow(endInd(1),endInd(2),2,k);
                endPixel = endPixel + [currFlowY, currFlowX];
            else
                k = (currPath(j+1)-1)*(currPath(j+1)-2)/2 + currPath(j);
                endInd = round(endPixel);
                currFlowX = -flow(endInd(1),endInd(2),1,k);
                currFlowY = -flow(endInd(1),endInd(2),2,k);
                endPixel = endPixel + [currFlowY, currFlowX];
            end
        end
        EndPixels(i,:) = endPixel;
    else
        EndPixels(i,:) = endPixel;
    end
end

DistEndPixels = sqrt(sum((EndPixels-repmat(expectedPixel,length(path),1)).^2,2));
[minDist, minInd] = min(DistEndPixels);
FinalPath{count} = path{minInd};
path = FindingPath(mov, DistMat, minInd, mean(DistMat(:)));
end

nNodes=1;
for i=1:length(FinalPath)
    for j=2:length(FinalPath{i})
        SingleFinalPath(nNodes) = FinalPath{i}(j);
        nNodes = nNodes+1;
    end
end
SingleFinalPath = [1 SingleFinalPath];
figure;
for i=1:length(SingleFinalPath)
    imshow(mov(:,:,:,SingleFinalPath(i)));
end

%Render Slow Motion Interpolation:
% PathMov = MotionInterp(mov, SingleFinalPath);
% for i=1:length(SingleFinalPath)-1
%     imshow(PathMov(:,:,:,i));
% end

boolContinue = input('Do you want to continue? (y/n): ','s');
if (boolContinue == 'n')
    break;
end
close all;
clear FinalPath;
clear SingleFinalPath;
end

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

function [DistMat1, ProbMat1] = DP_DistMat(mov)
%**************************************************************************
% Compute the Advanced Distance Matrix and Probability Matrix corresponding to 
% Dynamics Preserving:
[nRows nColumns nColors nFrames] = size(mov);
kernelWeights = [1 4 6 4 1]/16;
[DistMat ProbMat] = BasicDistMat(mov);
DistMat1 = DistMat;
ProbMat1 = ProbMat;
for j = 3:nFrames-2
    for i = 3:nFrames-2
        DistMat1(i,j) = kernelWeights*[DistMat(i-2,j-2) DistMat(i-1,j-1) DistMat(i,j) DistMat(i+1,j+1) DistMat(i+2,j+2)]';
    end
end

sigma1 = sum(DistMat1(:))/sum(sum(DistMat1~=0));
ProbMat1(1:end-1,:) = exp(-DistMat1(2:end,:)/sigma1);
sumVec1 = sum(ProbMat1(1:end-1,:),2);
ProbMat1(1:end-1,:) = ProbMat1(1:end-1,:)./repmat(sumVec1,1,nFrames);
end
%**************************************************************************

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

function flow = CreateFlow(mov)
[nRows nColumns nColors nFrames] = size(mov);
flow = zeros(nRows,nColumns,2,nFrames*(nFrames-1)/2);
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

fprintf('%d pieces of flows to be computed.\n',nFrames*(nFrames-1)/2);

k=1;

    for i = 2:20
        for j = 1:i-1
            im1 = mov(:,:,:,i);
            im2 = mov(:,:,:,j);
            [vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
            flow(:,:,1,k) = vx;
            flow(:,:,2,k) = vy;
            fprintf('%d is done.\n',k);
            k=k+1;
        end
    end
end

function PathMov = MotionInterp(mov, FinalPath)

alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;
para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

for i=1:length(FinalPath)-1
    [vxTemp vyTemp PathMov(:,:,:,i)] =...
    Coarse2FineTwoFrames(mov(:,:,:,FinalPath(i)),mov(:,:,:,FinalPath(i+1)),para);
end
end