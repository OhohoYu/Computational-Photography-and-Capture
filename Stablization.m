function newMov = Stablization( movOriginal,movGray, window_size)
[nRows, nColumns, nColors, nFrames] = size(movOriginal);

newMov = movOriginal;
Tx=zeros(1,nFrames);
Ty=zeros(1,nFrames);
Rtheta=zeros(1,nFrames);
Scaling=ones(1,nFrames);
for i=2:nFrames
    % Step 1. Read Frames from a Movie File
    % Since most of the static features are on the upper half of the frame
    % Here we only take part of frame to extract features:
    % NOTE: 'range' could be 1:nY
    imgA = movGray(:,:,i-1);
    FimgA = imgA(:,:);
    imgB = movGray(:,:,i);
    FimgB = imgB(:,:);
    
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
                                  pointsB, pointsA, 'similarity','MaxDistance',0.3);

    % Step 5. Transform Approximation and Smoothing
    % Extract scale and rotation part sub-matrix.
    H = tform.T;
    R = H(1:2,1:2);
    % Compute theta from mean of two possible arctangents
    theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
    % Compute scale from mean of two stable mean calculations
    scale = mean(R([1 4])/cos(theta));
    Scaling(i) = scale;
    % Translation remains the same:
    translation = H(3, 1:2);
    
    fprintf('processing the No.%d Frame\n',i);
    Tx(i)=translation(1);
    Ty(i)=translation(2);
    Rtheta(i)=theta;
end

TrajectX = cumsum(Tx);
TrajectY = cumsum(Ty);
TrajectTheta = cumsum(Rtheta);

%window_size = 5;
STX = tsmovavg(TrajectX,'s',window_size);
STY = tsmovavg(TrajectY,'s',window_size);
STT = tsmovavg(TrajectTheta,'s',window_size);

for i=1:window_size-1
    STX(i) = mean(TrajectX(i:i+window_size-1));
    STY(i) = mean(TrajectY(i:i+window_size-1));
    STT(i) = mean(TrajectTheta(i:i+window_size-1));
end

Tx = Tx-(STX-TrajectX);
Ty = Ty-(STY-TrajectY);
Rtheta = Rtheta-(STT-TrajectTheta);

for i=1:nFrames
    im = movOriginal(:,:,:,i);
    % Reconstitute new s-R-t transform:
    scale = Scaling(i);
    theta = Rtheta(i);
    translation = [Tx(i),Ty(i)];
    HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
            translation], [0 0 1]'];
    tformsRT = affine2d(HsRt);    
    imgBsRt = imwarp(im, tformsRT, 'OutputView', imref2d(size(im)));
    newMov(:,:,:,i) = imgBsRt;
    
end
end

