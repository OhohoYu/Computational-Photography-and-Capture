close all;
clear all;

% HorizontalCube = load_sequence_color('./Structured Light Data/cube_T1','00', 0, 19, 2, 'png');
% VerticalCube = load_sequence_color('./Structured Light Data/cube_T1','00', 20, 39, 2, 'png');
% 
% GrayHorizontalCube = zeros(2048,3072,20);
% GrayVerticalCube = zeros(2048,3072,20);
% 
% for i = 1:20
%     test1 = rgb2gray(HorizontalCube(:,:,:,i));
%     test2 = rgb2gray(VerticalCube(:,:,:,i));
%     GrayHorizontalCube(:,:,i) = test1;
%     GrayVerticalCube(:,:,i) = test2;
% end

GrayHorizontalCube = load('GrayHorizontalCube.mat');
GrayHorizontalCube = GrayHorizontalCube.GrayHorizontalCube;
GrayVerticalCube = load('GrayVerticalCube.mat');
GrayVerticalCube = GrayVerticalCube.GrayVerticalCube;

[nRows nColumns nFrames] = size(GrayVerticalCube);
nPairs =10;
threshold = 0.004;
for i=1:nPairs  
    temp1 =GrayHorizontalCube(:,:,2*i) - GrayHorizontalCube(:,:,2*i-1);
    temp1(find(abs(temp1)<threshold)) = 0;
    HorizontalDiff(:,:,i) = temp1;
   
    temp2= GrayVerticalCube(:,:,2*i) - GrayVerticalCube(:,:,2*i-1);
    temp2(find(abs(temp2)<threshold)) = 0;
    VerticalDiff(:,:,i) = temp2;
end

UV_Data = []; %projector
XY_Data = []; %camera
nPoints = 0;
for iRow = 1:nRows
    fprintf('Processing No.%d Row:\n',iRow);
    for iColumn = 1:1:nColumns
        if (find(HorizontalDiff(iRow,iColumn,:)==0))
            continue;
        end
        if (find(VerticalDiff(iRow,iColumn,:)==0))
            continue;
        end        
        currHoriDiff = reshape(HorizontalDiff(iRow,iColumn,:),1,nPairs);
        currVertDiff = reshape(VerticalDiff(iRow,iColumn,:),1,nPairs);
        currHoriBinaryCode = (currHoriDiff<0);
        currVertBinaryCode = (currVertDiff<0);
        u = bi2de(currVertBinaryCode,2,'right-msb')+1;
        v = bi2de(currHoriBinaryCode,2,'right-msb')+1;
        XY_Data = [XY_Data [iRow;iColumn]];
        UV_Data = [UV_Data [v;u]];
        nPoints = nPoints+1;
    end
end

CamIntrinsic = [   4786.25390625,       0.00000000,    1541.16491699;
                     0.00000000,    4789.81884766,    1036.94421387;
                     0.00000000,       0.00000000,       1.00000000];
ProjIntrinsic = [   3680.39404297,       0.00000000,     591.75494385;
                     0.00000000,    3672.32153320,     393.62173462;
                     0.00000000,       0.00000000,       1.00000000];
CamT = [-0.10157115; -0.10139455; 0.49625999];
ProjT = [-0.14915472; -0.06104425;  1.36014771];
CamR = [0.99998617,      -0.00475739,       0.00223672;
        0.00382316,       0.95016861,       0.31171292;
       -0.00360820,      -0.31170005,       0.95017368];
ProjR = [0.72119248,       0.44233182,      -0.53312665;
        -0.36164442,       0.89680630,       0.25485638;
         0.59084243,       0.00900178,       0.80673677];

% Reconstruction using my own codes implemented in lab2
%**************************************************************************
w = zeros(3,size(UV_Data,2));
ImPoints = zeros(size(UV_Data,2),4);
ImPoints(:,1) = UV_Data(2,:)'; % x-coord in projector
ImPoints(:,2) = UV_Data(1,:)'; % y-coord in projector
ImPoints(:,3) = XY_Data(2,:)'; % x-coord in camera
ImPoints(:,4) = XY_Data(1,:)'; % y-coord in camera

for i=1:size(UV_Data,2)
    A = [];
    b = [];
    for j=1:2
        switch j
            case 1
            R = ProjR;
            T = ProjT;
            JointsCam = ProjIntrinsic\[ImPoints(i,2*j-1),ImPoints(i,2*j),1]';
            case 2
            R = CamR;
            T = CamT;
            JointsCam = CamIntrinsic\[ImPoints(i,2*j-1),ImPoints(i,2*j),1]';
        end
        a1 = [R(3,1)*JointsCam(1)-R(1,1),R(3,2)*JointsCam(1)-R(1,2),R(3,3)*JointsCam(1)-R(1,3)];
        A = [A;a1];
        a2 = [R(3,1)*JointsCam(2)-R(2,1),R(3,2)*JointsCam(2)-R(2,2),R(3,3)*JointsCam(2)-R(2,3)];
        A = [A;a2];
        b1 = [T(1)-T(3)*JointsCam(1);T(2)-T(3)*JointsCam(2)];
        b = [b;b1];
    end
    w(:,i) = (A'*A)^(-1)*A'*b;
end
%**************************************************************************
P = w;

% Reconstruction Test Using Sample Solution of Lab2
%**************************************************************************
% p1 = [XY_Data(2,:);XY_Data(1,:);ones(1,nPoints)]; %Camera
% p2 = [UV_Data(2,:);UV_Data(1,:);ones(1,nPoints)]; %Projector
% 
% p1p = CamIntrinsic\p1; %equivalent to typing inv(intrinsic)*p1
% p2p = ProjIntrinsic\p2;
% 
% P = zeros(3,nPoints);
% 
% for i=1:nPoints
%     A = zeros(4,3);
%     b = zeros(4,1);
%     
%     % camera
%     A(1,:) = [CamR(3,1)*p1p(1,i)-CamR(1,1), CamR(3,2)*p1p(1,i)-CamR(1,2), CamR(3,3)*p1p(1,i)-CamR(1,3)];
%     A(2,:) = [CamR(3,1)*p1p(2,i)-CamR(2,1), CamR(3,2)*p1p(2,i)-CamR(2,2), CamR(3,3)*p1p(2,i)-CamR(2,3)];
%     
%     b(1:2) = [CamT(1)-CamT(3)*p1p(1,i); CamT(2)-CamT(3)*p1p(2,i)];
%     
%     % projector
%     A(3,:) = [ProjR(3,1)*p2p(1,i)-ProjR(1,1), ProjR(3,2)*p2p(1,i)-ProjR(1,2), ProjR(3,3)*p2p(1,i)-ProjR(1,3)];
%     A(4,:) = [ProjR(3,1)*p2p(2,i)-ProjR(2,1), ProjR(3,2)*p2p(2,i)-ProjR(2,2), ProjR(3,3)*p2p(2,i)-ProjR(2,3)];
%     
%     b(3:4) = [ProjT(1)-ProjT(3)*p2p(1,i); ProjT(2)-ProjT(3)*p2p(2,i)];
%     
%     %compute the least squares solution
%     w = A\b; %equivalent to w=inv(A'*A)*A'*b
%     
%     %store the computed point in an array
%     P(:,i) = w;
%     
% end
%**************************************************************************

scatter3(P(1,:),P(2,:),P(3,:));


maxP = max(P(3,:));
minP = min(P(3,:));
P(3,:) = (maxP-P(3,:))./(maxP-minP);

depthMap = zeros(nRows,nColumns);

for i=1:nPoints-1
    depthMap(XY_Data(1,i),XY_Data(2,i))=P(3,i);
    depthMap(XY_Data(1,i+1),XY_Data(2,i+1))=P(3,i+1);
    if (XY_Data(1,i)==XY_Data(1,i+1))
        nPixels = XY_Data(2,i+1)-XY_Data(2,i)-1;
        for j=1:nPixels
            depthMap(XY_Data(1,i),XY_Data(2,i)+j)=P(3,i)+...
                j*(P(3,i+1)-P(3,i))/(nPixels+1);
        end
    end
end
figure;
imshow(depthMap);

% Create a ground truth for correction
%**************************************************************************
% im1 = imread('./Structured Light Data/cube_T1/0039.png');
% im2 = imread('./Structured Light Data/cube_T1/0038.png');
% im1 = im2double(im1);
% im2 = im2double(im2);
% im1 = rgb2gray(im1);
% im2 = rgb2gray(im2);
% GroundTruth = im1+im2;
% GroundTruth(GroundTruth>0.5) = 1;
% GroundTruth(GroundTruth<=0.5) = 0;
%**************************************************************************
GroundTruth = load('GroundTruth.mat');
GroundTruth = GroundTruth.GroundTruth;

testSet = depthMap;
testSet(testSet>0) = 1;
testSet(testSet<=0) = 0;

correctedDepthMap = depthMap;
for i=1:nRows
    for j=1:nColumns
        if (GroundTruth(i,j)==1 && depthMap(i,j)==0)
            mask = CreateMask(testSet,i,j,15);
            maskIm = CreateMask(depthMap,i,j,15);
            maskIm(isnan(maskIm))=0;
            correctedDepthMap(i,j) = sum(sum(maskIm.*(mask)))/sum(mask(:));
        end
    end
end
figure;
imshow(correctedDepthMap);

function mask = CreateMask(Matrix,i,j,size)
mask = zeros(size,size);
for n = 1:size 
    for m = 1:size
       mask(n,m) = Matrix(i-(size-1)/2+(n-1),j-(size-1)/2+(m-1));
    end
end
end