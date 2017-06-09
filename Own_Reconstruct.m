close all;
clear all;

GrayHorizontalCube = load_sequence_color('./Data','OWN_00', 00, 19, 2, 'JPG');
GrayVerticalCube = load_sequence_color('./Data','OWN_00', 20, 39, 2, 'JPG');

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

CamIntrinsic = [3918.50023416561,0,1125.56714177458;0,3897.31634257115,826.85724068224;0,0,1];
ProjIntrinsic =[12870.40541,0,1190.08440;0,13142.76551,833.29423;0,0,1];
CamT = [ -193.098592 	; -170.863088 	; 3406.097517 ];
ProjT = [ -257.789727 ;	 -245.397909 	; 5146.529704 ];
CamR = [ 0.023505 	 0.981246 	 0.292285;
                              0.998120 	 -0.074727 	 0.170083;
                              0.189785 	 0.281167 	 -0.979735 ];
ProjR = [ -0.015695 	 0.998929 	 -0.217252;
                               0.966981 	 -0.026053 	 -0.180652;
                               -0.180718 	 -0.224441 	 -0.968050 ];


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
