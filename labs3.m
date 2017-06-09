function output = labs3%(path, prefix, first, last, digits, suffix)

%
% Read a sequence of images and correct the film defects. This is the file 
% you have to fill for the coursework. Do not change the function 
% declaration, keep this skeleton. You are advised to create subfunctions.
% 
% Arguments:
%
% path: path of the files
% prefix: prefix of the filename
% first: first frame
% last: last frame
% digits: number of digits of the frame number
% suffix: suffix of the filename
%
% This should generate corrected images named [path]/corrected_[prefix][number].png
%
% Example:
%
% mov = labs3('../images','myimage', 0, 10, 4, 'png')
%   -> that will load and correct images from '../images/myimage0000.png' to '../images/myimage0010.png'
%   -> and export '../images/corrected_myimage0000.png' to '../images/corrected_myimage0010.png'
%

% mov3 = load_sequence('../lab3/corrected3','myimage', 1, 657, 3, 'png');
% mov4 = load_sequence('../lab3/corrected4','myimage', 1, 657, 3, 'png');

%Load Frames
mov = load_sequence('../lab3/footage','footage_', 1, 657, 3, 'png');
mov = im2double(mov);
[nY,nX,nFrames] = size(mov);

%==========================================================================
%==========================================================================

%Part 1
%Detection of Scene Cuts

%Compute Average Greyscale and Approx Difference
%Find the Number of the Frames that have larger Difference,
%which denote the scene cuts
aveGrey = mean(mean(mov,1),2);
aveGrey = reshape(aveGrey,1,nFrames);
DeriAveGrey = diff(aveGrey);
diffThresh = 0.08;
%Record the number of scene cuts and index
nCuts = 0;
cutsIndex = [];
for i=1:nFrames-1
    if abs(DeriAveGrey(i)) > diffThresh
        nCuts = nCuts+1;
        cutsIndex = [cutsIndex i];
    end
end
nScenes = nCuts+1;
sceneIndex = zeros(2,nCuts+1);
sceneIndex(:,1) = [1;cutsIndex(1)];
for i=2:nCuts
    sceneIndex(:,i) = [cutsIndex(i-1)+1;cutsIndex(i)];
end
sceneIndex(:,end) = [cutsIndex(nCuts)+1;nFrames];

%Output movie1.avi for part 1
% mov1 = mov;
% v1 = VideoWriter('movie1.avi');
% tDuration = 10;
% v1.FrameRate = nFrames/tDuration;
% open(v1);
% for i=1:nFrames
%     A = mov1(:,:,i);
%     for j=1:nCuts
%         if i >= cutsIndex(j)+1 && i <= cutsIndex(j)+20
%            A = insertText(A, [20, 20], 'change');
%         end
%     end
%     writeVideo(v1,A);
% end
% close(v1);

%==========================================================================
%==========================================================================

%Part 2
%Correction of Global Flicker:
Flickers = findFlickers(DeriAveGrey,cutsIndex);
%Add some missed flickers, which last more than 6 frames
Flickers = [Flickers [82;108] [178;199]];
%Compute corresponding avegrage intensity of flickers' neighbors:
FlickerIntst = NearAveInt(Flickers,cutsIndex,aveGrey);
%Output movie2.avi and all corrected frames of part 2
%Correct the flickering frames
%Use neighbor's average intensity to subsitute flickering frame's
nFlickers = size(Flickers,2);
% mov2 = mov;
% for i=1:nFlickers
%     startFrame = Flickers(1,i);
%     endFrame = Flickers(2,i);
%     for j=startFrame:endFrame
%         A = mov2(:,:,j);
%         mov2(:,:,j) = A*FlickerIntst(i)/(sum(A(:))/(nY*nX));
%     end
% end
% v2 = VideoWriter('movie2.avi');
% tDuration = 10;
% v2.FrameRate = nFrames/tDuration;
% open(v2);
% for i=1:nFrames
%     A = mov2(:,:,i);
%     writeVideo(v2,A);
% end
% close(v2);
% % save the output frames
% save_sequence(mov2, './corrected2', 'myimage', 1, 3);


%==========================================================================
%==========================================================================

%Part 3
%Correction of Blotches:
%Load Frames
% mov3 = load_sequence('../lab3/corrected2','myimage', 1, 657, 3, 'png');
% mov3 = im2double(mov3);
% [nY,nX,nFrames] = size(mov3);
% tempMov = mov3;
% 
% for i = 1:nScenes
%     for j = sceneIndex(1,i)+1:sceneIndex(2,i)-1
%         fprintf('processing the No.%d Frame\n',j);
%         preIm = tempMov(:,:,j-1);
%         currIm = tempMov(:,:,j);
%         nextIm = tempMov(:,:,j+1);
%         mov3(:,:,j) = Blotches(preIm,currIm,nextIm);
%     end
% end
% 
% v3 = VideoWriter('movie3.avi');
% tDuration = 10;
% v3.FrameRate = nFrames/tDuration;
% open(v3);
% for i=1:nFrames
%     A = mov3(:,:,i);
%     writeVideo(v3,A);
% end
% close(v3);
% % save the output frames
% save_sequence(mov3, './corrected3', 'myimage', 1, 3);

%==========================================================================
%==========================================================================

%Part 4
%Correction of vertical artefacts:
%Load Frames
mov4 = load_sequence('../lab3/corrected3','myimage', 1, 657, 3, 'png');
mov4 = im2double(mov4);
[nY,nX,nFrames] = size(mov4);

v4 = VideoWriter('movie4.avi');
tDuration = 10;
v4.FrameRate = nFrames/tDuration;
open(v4);

for i = 1:nFrames
    A = mov4(:,:,i);
    if i>=sceneIndex(1,end) && i<=sceneIndex(2,end)
        A = verStripes(A);
        writeVideo(v4,A);
        mov4(:,:,i) = A;
    else
        writeVideo(v4,A);
        mov4(:,:,i) = A;
    end
end
close(v4);
% save the output frames
save_sequence(mov4, './corrected4', 'myimage', 1, 3);
%==========================================================================
%==========================================================================

%Part 5
%Correction of camera shake :
%Load Frames
% mov5 = load_sequence('../lab3/corrected4','myimage', 1, 657, 3, 'png');
% mov5 = im2double(mov5);
% 
% scene1 = mov5(:,:,sceneIndex(1,1):sceneIndex(2,1));
% newScene1 = VideoStabil(scene1,300,'similarity',0.5);
% 
% scene2 = mov5(:,:,sceneIndex(1,2):sceneIndex(2,2));
% newScene2 = VideoStabil(scene2,300,'similarity',0.5);
% 
% scene3 = mov5(:,:,sceneIndex(1,3):sceneIndex(2,3));
% newScene3 = VideoStabil(scene3,300,'similarity',0.5);
% 
% %average frame:
% %imshow(sum(mov,3)/nFrames);
% 
% mov5 = cat(3,newScene1,newScene2,newScene3);
% 
% v5 = VideoWriter('movieTest.avi');
% tDuration = 10;
% v5.FrameRate = nFrames/tDuration;
% open(v5);
% for i=1:161
%     A = mov5(:,:,i);
%     writeVideo(v5,A);
% end
% close(v5);
% save_sequence(mov5, './corrected5', 'myimage', 1, 3);




%Subfunctions:
%==========================================================================
%Flickers Detection
function Flickers = findFlickers(DeriAveGrey,cutsIndex)
% plot(DeriAveGrey);
% Manually set Threshold=0.02, to find the suspected flickering frames,
% whose average intensity difference is relatively large
Flickers = [];
flickersIndex = [];
for i=1:size(DeriAveGrey,2)
    if abs(DeriAveGrey(i)) > 0.02
        flickersIndex = [flickersIndex i];
    end
end

%Find the starting and ending frames for each flickers
%Here we define some conditions for searching:
%Flicker lasts within 6 frames.
%The differences of the starting and ending frame have different signs
for i=1:length(flickersIndex)-1
    for j=length(flickersIndex):-1:i+1
        if flickersIndex(j)-flickersIndex(i)<=6 && flickersIndex(j)-flickersIndex(i)>1
            if DeriAveGrey(flickersIndex(j))*DeriAveGrey(flickersIndex(i)) < 0
                Flickers = [Flickers [flickersIndex(i)+1;flickersIndex(j)]];
                break;
            end
        end
    end
end

%Eliminate the innocent "flicker", which are actually in two scenes.
tempFlickers = Flickers;
for i=1:size(Flickers,2)
    for j=1:length(cutsIndex)
        if cutsIndex(j)>=Flickers(1,i) && cutsIndex(j)<=Flickers(2,i)
            tempFlickers(:,i)=[];
        end
    end
end
Flickers = tempFlickers;

%==========================================================================
%Flickers Intensity Correction
function FlickerIntst = NearAveInt(Flickers,cutsIndex,aveGrey)
nFlickers = size(Flickers,2);
FlickerIntst = zeros(1,nFlickers);
%Compute the average intensity of flickers' neighbors
for i=1:size(Flickers,2)
    for j=1:length(cutsIndex)
        if Flickers(1,i)<=(cutsIndex(j)+5) && Flickers(1,i)>=(cutsIndex(j)-5)
           FlickerIntst(i) = mean(aveGrey(Flickers(2,i)+1:Flickers(2,i)+10)); 
           break;
        elseif Flickers(2,i)<=(cutsIndex(j)+5) && Flickers(2,i)>=(cutsIndex(j)-5)
                FlickerIntst(i) = mean(aveGrey(Flickers(1,i)-10:Flickers(1,i)-1));
           break;
        else
            FlickerIntst(i) = (sum(aveGrey(Flickers(1,i)-5:Flickers(1,i)-1))...
                               +sum(aveGrey(Flickers(2,i)+1:Flickers(2,i)+5)))/10;
           break;
        end
    end
end
%==========================================================================
%Detection and Correction of Blotches:
function NewCurr = Blotches(imPre,imCurr,imNext)
%Blotches Detection:
NewCurr = imCurr;
imDiff1 = abs(imPre-imCurr);
imDiff2 = abs(imCurr-imNext);
imDiff1(imDiff1>0.1)=1;
imDiff1(imDiff1<0.1)=0;
imDiff2(imDiff2>0.1)=1;
imDiff2(imDiff2<0.1)=0;
imBlo=imDiff1.*imDiff2;
%Blotches Correction: Create a 3-by-3 mask
[imY,imX] = size(imBlo);
for i = 2:imY-1
    for j = 2:imX-1
        if imBlo(i,j) == 1
            maskBlo = [imBlo(i-1,j-1) imBlo(i-1,j) imBlo(i-1,j+1);...
                       imBlo(i,j-1) imBlo(i,j) imBlo(i,j+1);...
                       imBlo(i+1,j-1) imBlo(i+1,j) imBlo(i+1,j+1)];
            maskIm = [NewCurr(i-1,j-1) NewCurr(i-1,j) NewCurr(i-1,j+1);...
                       NewCurr(i,j-1) NewCurr(i,j) NewCurr(i,j+1);...
                       NewCurr(i+1,j-1) NewCurr(i+1,j) NewCurr(i+1,j+1)];
            NewCurr(i,j) = sum(sum(maskIm.*(~maskBlo)))/sum(~maskBlo(:));
            imBlo(i,j) = 0;
        end
    end
end

%==========================================================================
%Vertical Stripes Removal:
function NewCurr = verStripes(imCurr)
%Remove vertical stripes using median filter:
[nY,nX]=size(imCurr);
temp = reshape(imCurr',1,nY*nX);
temp = medfilt1(temp);
temp = reshape(temp,nX,nY);
NewCurr = temp';

%==========================================================================
%:
function newMov = VideoStabil(mov,range,TransformType,maxDist)
%Video Stabilization Using Point Feature Matching
[nY nX nFrames] = size(mov);
newMov = mov;
for i=2:nFrames 
    % Step 1. Read Frames from a Movie File
    % Since most of the static features are on the upper half of the frame
    % Here we only take part of frame to extract features:
    % NOTE: 'range' could be 1:nY
    imgA = newMov(:,:,i-1);
    FimgA = imgA(1:range,:);
    imgB = newMov(:,:,i);
    FimgB = imgB(1:range,:);
    
    % Step 2. Collect SURF Points from Each Frame
    pointsA = detectSURFFeatures(FimgA);
    pointsB = detectSURFFeatures(FimgB);
    
    % Step 3. Select Correspondences Between Points
    % Extract FREAK descriptors for the corners
    [featuresA, pointsA] = extractFeatures(imgA, pointsA);
    [featuresB, pointsB] = extractFeatures(imgB, pointsB);

    indexPairs = matchFeatures(featuresA, featuresB);
    pointsA = pointsA(indexPairs(:, 1), :);
    pointsB = pointsB(indexPairs(:, 2), :);

    % Step 4. Estimating Transform from Noisy Correspondences
    % TransformType could be 'similarity','affine', 'projective'
    [tform, pointsBm, pointsAm] = estimateGeometricTransform(...
                                  pointsB, pointsA, TransformType,'MaxDistance',maxDist);

    % Step 5. Transform Approximation and Smoothing
    % Extract scale and rotation part sub-matrix.
    H = tform.T;
    R = H(1:2,1:2);
    % Compute theta from mean of two possible arctangents
    theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
    % Compute scale from mean of two stable mean calculations
    scale = mean(R([1 4])/cos(theta));
    % Translation remains the same:
    translation = H(3, 1:2);
    % Reconstitute new s-R-t transform:
    HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
            translation], [0 0 1]'];
    tformsRT = affine2d(HsRt);    
    imgBsRt = imwarp(imgB, tformsRT, 'OutputView', imref2d(size(imgB)));
    
    fprintf('processing the No.%d Frame\n',i);
    
    % Add colors to the 'empty'(black) area of the output frame
    patch = sum(mov(:,:,1:i),3)/i;
    imgBsRt(find(imgBsRt==0)) = patch(find(imgBsRt==0));
    
    newMov(:,:,i) = imgBsRt;
    
end