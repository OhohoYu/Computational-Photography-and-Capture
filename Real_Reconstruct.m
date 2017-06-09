close all;
clear all;

GrayHorizontalCube = load_sequence_color('./Structured Light Data/real_crayon_dalek','IMG_94', 18, 37, 2, 'JPG');
GrayVerticalCube = load_sequence_color('./Structured Light Data/real_crayon_dalek','IMG_94', 38, 57, 2, 'JPG');

% GrayHorizontalCube = zeros(2048,3072,20);
% GrayVerticalCube = zeros(2048,3072,20);
% 
% for i = 1:20
%     test1 = rgb2gray(HorizontalCube(:,:,:,i));
%     test2 = rgb2gray(VerticalCube(:,:,:,i));
%     GrayHorizontalCube(:,:,i) = test1;
%     GrayVerticalCube(:,:,i) = test2;
% end

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
    for iColumn = 1:4:nColumns
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

CamIntrinsic = [3894.60930532672,0,2495.75038547156;0,3980.10145847271,2256.78576061699;0,0,1];
ProjIntrinsic =[1159.35664942873,0,511.500000000000;0,957.521815550119,383.500000000000;0,0,1];
CamT = [ -510.804704; -597.597144; 2026.343659 ];
ProjT = [ -103.029020 ; -32.393755 ; 586.754791 ];
CamR = [ -0.089925 	 0.990768 	 -0.101453;
        0.974874 	 0.066719 	 -0.212531;
         -0.203800 	 -0.118016 	 -0.971873 ];
ProjR = [ 0.263801 	 0.912316 	 0.302739;
         0.809270 	 -0.062374 	 -0.573254;
        -0.490132 	0.418123       -0.732642 ];


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
