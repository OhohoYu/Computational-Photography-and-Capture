clear all;
close all;
mov = load_sequence_color('./Drones','testFootage', 0, 833, 3, 'jpg');

Ind = load('Ind.mat');
Ind = Ind.Ind;

Gray = load('movGray.mat');
Gray = Gray.movGray;
[nY nX nFrames] = size(Gray);

movGray = zeros(nY,nX,length(Ind));
movOriginal = zeros(nY,nX,3,length(Ind));

for i=1:length(Ind)
    movGray(:,:,i) = Gray(:,:,Ind(i));
    movOriginal(:,:,:,i) = mov(:,:,:,Ind(i));
end

[nRows nColumns nFrames] = size(movGray);
% newMov = movOriginal;

mov1 = movOriginal(:,:,:,1:161);
gray1 = movGray(:,:,1:161);
new1 = Stablization(mov1,gray1,20);

mov2 = movOriginal(:,:,:,191:309);
gray2 = movGray(:,:,191:309);
new2 = Stablization(mov2,gray2,10);

mov3 = movOriginal(:,:,:,361:384);
gray3 = movGray(:,:,361:384);
new3 = Stablization(mov3,gray3,5);

mov4 = movOriginal(:,:,:,406:447);
gray4 = movGray(:,:,406:447);
new4 = Stablization(mov4,gray4,5);

mov5 = movOriginal(:,:,:,457:494);
gray5 = movGray(:,:,457:494);
new5 = Stablization(mov5,gray5,5);

mov6 = movOriginal(:,:,:,505:end);
gray6 = movGray(:,:,505:end);
new6 = Stablization(mov6,gray6,5);

newMov = cat(4,new1,movOriginal(:,:,:,162:190));
newMov = cat(4,newMov,new2);
newMov = cat(4,newMov,movOriginal(:,:,:,310:360));
newMov = cat(4,newMov,new3);
newMov = cat(4,newMov,movOriginal(:,:,:,385:405));
newMov = cat(4,newMov,new4);
newMov = cat(4,newMov,movOriginal(:,:,:,448:456));
newMov = cat(4,newMov,new5);
newMov = cat(4,newMov,movOriginal(:,:,:,495:504));
newMov = cat(4,newMov,new6);

% newMov(:,:,:,1:161) = new1;
% newMov(:,:,:,191:309) = new2;
% newMov(:,:,:,361:384) = new3;
% newMov(:,:,:,406:447) = new4;
% newMov(:,:,:,457:494) = new5;
% newMov(:,:,:,505:end) = new6;

% [nRows, nColumns, nColors, nFrames] = size(movOriginal);
% 
% newMov = movOriginal;
% Tx=zeros(1,nFrames);
% Ty=zeros(1,nFrames);
% Rtheta=zeros(1,nFrames);
% Scaling=ones(1,nFrames);
% for i=2:nFrames
%     % Step 1. Read Frames from a Movie File
%     % Since most of the static features are on the upper half of the frame
%     % Here we only take part of frame to extract features:
%     % NOTE: 'range' could be 1:nY
%     imgA = movGray(:,:,i-1);
%     FimgA = imgA(:,:);
%     imgB = movGray(:,:,i);
%     OriginalB = movOriginal(:,:,:,i);
%     FimgB = imgB(:,:);
%     
%     % Step 2. Collect SURF Points from Each Frame
%     pointsA = detectSURFFeatures(FimgA);
%     pointsB = detectSURFFeatures(FimgB);
%     
%     % Step 3. Select Correspondences Between Points
%     % Extract FREAK descriptors for the corners
%     [featuresA, pointsA] = extractFeatures(imgA, pointsA);
%     [featuresB, pointsB] = extractFeatures(imgB, pointsB);
% 
%     indexPairs = matchFeatures(featuresA, featuresB);
%     pointsA = pointsA(indexPairs(:, 1), :);
%     pointsB = pointsB(indexPairs(:, 2), :);
% 
%     % Step 4. Estimating Transform from Noisy Correspondences
%     % TransformType could be 'similarity','affine', 'projective'
%     [tform, pointsBm, pointsAm] = estimateGeometricTransform(...
%                                   pointsB, pointsA, 'similarity','MaxDistance',0.3);
% 
%     % Step 5. Transform Approximation and Smoothing
%     % Extract scale and rotation part sub-matrix.
%     H = tform.T;
%     R = H(1:2,1:2);
%     % Compute theta from mean of two possible arctangents
%     theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
%     % Compute scale from mean of two stable mean calculations
%     scale = mean(R([1 4])/cos(theta));
%     Scaling(i) = scale;
%     % Translation remains the same:
%     translation = H(3, 1:2);
%     
%     fprintf('processing the No.%d Frame\n',i);
%     Tx(i)=translation(1);
%     Ty(i)=translation(2);
%     Rtheta(i)=theta;
% end
% 
% TrajectX = cumsum(Tx);
% TrajectY = cumsum(Ty);
% TrajectTheta = cumsum(Rtheta);
% 
% window_size = 5;
% STX = tsmovavg(TrajectX,'s',window_size);
% STY = tsmovavg(TrajectY,'s',window_size);
% STT = tsmovavg(TrajectTheta,'s',window_size);
% 
% for i=1:window_size-1
%     STX(i)=0;
%     STY(i)=0;
%     STT(i)=0;
% end
% 
% Tx = Tx-(STX-TrajectX);
% Ty = Ty-(STY-TrajectY);
% Rtheta = Rtheta-(STT-TrajectTheta);
% 
% for i=1:nFrames
%     im = movOriginal(:,:,:,i);
%     % Reconstitute new s-R-t transform:
%     scale = Scaling(i);
%     theta = Rtheta(i);
%     translation = [Tx(i),Ty(i)];
%     HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
%             translation], [0 0 1]'];
%     tformsRT = affine2d(HsRt);    
%     imgBsRt = imwarp(im, tformsRT, 'OutputView', imref2d(size(im)));
%     newMov(:,:,:,i) = imgBsRt;
%     
% end

v1 = VideoWriter('movie1.avi');
v1.FrameRate = 24;
open(v1);
count=0;
for k=1:nFrames
    A1=rgb2gray(newMov(:,:,:,k));
    if (sum(~A1(:))/(nRows*nColumns))>0.3%||abs(Rtheta(k))>0.05
        A = movOriginal(:,:,:,k);
        count = count+1;
    else
        A = newMov(:,:,:,k);
    end
    A = A(:,:,:,:);
    writeVideo(v1,A);
    frame = newMov(:,:,:,k);
    if k<10
    imwrite(frame,strcat('Footage00',num2str(k),'.jpg'),'jpg');
    end
    if k<100&&k>=10
    imwrite(frame,strcat('Footage0',num2str(k),'.jpg'),'jpg');
    end
    if k<1000&&k>=100
    imwrite(frame,strcat('Footage',num2str(k),'.jpg'),'jpg');
    end
end
close(v1);
