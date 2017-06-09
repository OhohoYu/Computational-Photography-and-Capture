clear all;
close all;
movOriginal = load_sequence_color('./Drones','testFootage', 0, 833, 3, 'jpg');
[nRows, nColumns, nColors, nFrames] = size(movOriginal);

movGray = load('movGray.mat');
movGray = movGray.movGray;

avgMov = mean(mean(movGray,1),2);
avgMov = reshape(avgMov,1,nFrames);
diffMov = diff(avgMov);
CutPts = [1];
for i=1:length(diffMov)
    if abs(diffMov(i))>0.01
        CutPts = [CutPts i];
    end
end
CutPts = [CutPts nFrames];


Samples = cell((length(CutPts)-1),1);
fprintf('Segments: %d\n',(length(CutPts)-1));
count = 1;
for i=1:length(CutPts)-1
    startFrame = CutPts(i);
    endFrame = CutPts(i+1);
%     if (endFrame-startFrame)<=10
%         Samples{count} = startFrame:endFrame;
%         count = count+1;
%         fprintf('processing %d\n',count);
%         continue;
%     end    
    sceneMov = movOriginal(:,:,:,startFrame:endFrame);
    sceneGrayMov = movGray(:,:,startFrame:endFrame);
    [DistMat, ProbMat] = BasicDistMat(sceneMov);
    Threshold = mean(DistMat(:))-0.2;
    %Threshold = 0.1;
    Paths = FindingPath(sceneMov, DistMat, 1, Threshold);
    path = sort(Paths{length(Paths)});
    Samples{count} = path-1+startFrame;
    count = count+1;
    fprintf('processing %d\n',count);
end

Ind = [1];
for i=1:length(CutPts)-1
    currSamples = Samples{i}; 
    Ind = [Ind currSamples(2:end)];
end









