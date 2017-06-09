close all;
clear all;

% GrayHorizontalCube = load_sequence_color('./Structured Light Data/red_T1','00', 0, 19, 2, 'png');
% GrayVerticalCube = load_sequence_color('./Structured Light Data/red_T1','00', 20, 39, 2, 'png');
% 
% % GrayHorizontalCube = zeros(2048,3072,20);
% % GrayVerticalCube = zeros(2048,3072,20);
% % 
% % for i = 1:20
% %     test1 = rgb2gray(HorizontalCube(:,:,:,i));
% %     test2 = rgb2gray(VerticalCube(:,:,:,i));
% %     GrayHorizontalCube(:,:,i) = test1;
% %     GrayVerticalCube(:,:,i) = test2;
% % end
% 
% [nRows nColumns nFrames] = size(GrayVerticalCube);
% nPairs =10;
% threshold = 0.004;
% for i=1:nPairs  
%     temp1 =GrayHorizontalCube(:,:,2*i) - GrayHorizontalCube(:,:,2*i-1);
%     temp1(find(abs(temp1)<threshold)) = 0;
%     HorizontalDiff(:,:,i) = temp1;
%    
%     temp2= GrayVerticalCube(:,:,2*i) - GrayVerticalCube(:,:,2*i-1);
%     temp2(find(abs(temp2)<threshold)) = 0;
%     VerticalDiff(:,:,i) = temp2;
% end
% 
% UV_Data = []; %projector
% XY_Data = []; %camera
% nPoints = 0;
% for iRow = 1:nRows
%     fprintf('Processing No.%d Row:\n',iRow);
%     for iColumn = 1:5:nColumns
%         if (find(HorizontalDiff(iRow,iColumn,:)==0))
%             continue;
%         end
%         if (find(VerticalDiff(iRow,iColumn,:)==0))
%             continue;
%         end        
%         currHoriDiff = reshape(HorizontalDiff(iRow,iColumn,:),1,nPairs);
%         currVertDiff = reshape(VerticalDiff(iRow,iColumn,:),1,nPairs);
%         currHoriBinaryCode = (currHoriDiff<0);
%         currVertBinaryCode = (currVertDiff<0);
%         u = bi2de(currVertBinaryCode,2,'right-msb')+1;
%         v = bi2de(currHoriBinaryCode,2,'right-msb')+1;
%         XY_Data = [XY_Data [iRow;iColumn]];
%         UV_Data = [UV_Data [v;u]];
%         nPoints = nPoints+1;
%     end
% end
% 
% CamIntrinsic = [   4786.25390625,       0.00000000,    1541.16491699;
%                      0.00000000,    4789.81884766,    1036.94421387;
%                      0.00000000,       0.00000000,       1.00000000];
% ProjIntrinsic = [   3680.39404297,       0.00000000,     591.75494385;
%                      0.00000000,    3672.32153320,     393.62173462;
%                      0.00000000,       0.00000000,       1.00000000];
% CamT = [-0.10157115; -0.10139455; 0.49625999];
% ProjT = [-0.14915472; -0.06104425;  1.36014771];
% CamR = [0.99998617,      -0.00475739,       0.00223672;
%         0.00382316,       0.95016861,       0.31171292;
%        -0.00360820,      -0.31170005,       0.95017368];
% ProjR = [0.72119248,       0.44233182,      -0.53312665;
%         -0.36164442,       0.89680630,       0.25485638;
%          0.59084243,       0.00900178,       0.80673677];
% 
% % Reconstruction using my own codes implemented in lab2
% %**************************************************************************
% w = zeros(3,size(UV_Data,2));
% ImPoints = zeros(size(UV_Data,2),4);
% ImPoints(:,1) = UV_Data(2,:)'; % x-coord in projector
% ImPoints(:,2) = UV_Data(1,:)'; % y-coord in projector
% ImPoints(:,3) = XY_Data(2,:)'; % x-coord in camera
% ImPoints(:,4) = XY_Data(1,:)'; % y-coord in camera
% 
% for i=1:size(UV_Data,2)
%     A = [];
%     b = [];
%     for j=1:2
%         switch j
%             case 1
%             R = ProjR;
%             T = ProjT;
%             JointsCam = ProjIntrinsic\[ImPoints(i,2*j-1),ImPoints(i,2*j),1]';
%             case 2
%             R = CamR;
%             T = CamT;
%             JointsCam = CamIntrinsic\[ImPoints(i,2*j-1),ImPoints(i,2*j),1]';
%         end
%         a1 = [R(3,1)*JointsCam(1)-R(1,1),R(3,2)*JointsCam(1)-R(1,2),R(3,3)*JointsCam(1)-R(1,3)];
%         A = [A;a1];
%         a2 = [R(3,1)*JointsCam(2)-R(2,1),R(3,2)*JointsCam(2)-R(2,2),R(3,3)*JointsCam(2)-R(2,3)];
%         A = [A;a2];
%         b1 = [T(1)-T(3)*JointsCam(1);T(2)-T(3)*JointsCam(2)];
%         b = [b;b1];
%     end
%     w(:,i) = (A'*A)^(-1)*A'*b;
% end
% %**************************************************************************
% P = w;

UV_Data = load('monkey_UV.mat');
UV_Data = UV_Data.UV_Data;
XY_Data = load('monkey_XY.mat');
XY_Data = XY_Data.XY_Data;
P = load('monkey_P.mat');
P = P.P;

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
     
Camera = -CamR'*CamT;
CamExtrinsic = [[CamR CamT];0 0 0 1];
ProjExtrinsic = [[ProjR ProjT];0 0 0 1];
xImCart = projectiveCamera(CamIntrinsic,CamExtrinsic,P);
xImCart = round(xImCart);
Depth = P-repmat(Camera,1,(size(P,2)));
Depth = sqrt(Depth(1,:).^2+Depth(2,:).^2+Depth(3,:).^2);
minD = min(Depth);
maxD = max(Depth);
NormDepth = (maxD-Depth)./(maxD-minD);

nRows = 3048;
nColumns = 3072;
depthMap2 = zeros(nRows,nColumns);
% for i=1:size(P,2)-1
%     if (XY_Data(1,i)>=1&&XY_Data(1,i)<=2048&&XY_Data(1,i+1)>=1&&XY_Data(1,i+1)<=2048)
%     if (XY_Data(2,i)>=1&&XY_Data(2,i)<=3072&&XY_Data(2,i+1)>=1&&XY_Data(2,i+1)<=3072)
%     depthMap2(XY_Data(1,i),XY_Data(2,i))=NormDepth(i);
%     depthMap2(XY_Data(1,i+1),XY_Data(2,i+1))=NormDepth(i+1);
%     if (XY_Data(1,i)==XY_Data(1,i+1))
%         nPixels = XY_Data(2,i+1)-XY_Data(2,i)-1;
%         for j=1:nPixels
%             depthMap2(XY_Data(1,i),XY_Data(2,i)+j)=NormDepth(i)+...
%                 j*(NormDepth(i+1)-NormDepth(i))/(nPixels+1);
%         end
%     end
%     end
%     end
% end
xImCart = [xImCart(2,:);xImCart(1,:);linspace(1,size(xImCart,2),size(xImCart,2))];
xImCart = sortrows(xImCart')';
for i=1:size(P,2)-1
    if (xImCart(1,i)>=1&&xImCart(1,i)<=nRows&&xImCart(1,i+1)>=1&&xImCart(1,i+1)<=nRows)
    if (xImCart(2,i)>=1&&xImCart(2,i)<=nColumns&&xImCart(2,i+1)>=1&&xImCart(2,i+1)<=nColumns)
    depthMap2(xImCart(1,i),xImCart(2,i))=NormDepth(xImCart(3,i));
    depthMap2(xImCart(1,i+1),xImCart(2,i+1))=NormDepth(xImCart(3,i+1));
    if (xImCart(1,i)==xImCart(1,i+1))
        if (abs(xImCart(2,i+1)-xImCart(2,i))<=20)
        nPixels = xImCart(2,i+1)-xImCart(2,i)-1;
        for j=1:nPixels
            depthMap2(xImCart(1,i),xImCart(2,i)+j)=NormDepth(xImCart(3,i))+...
                j*(NormDepth(xImCart(3,i+1))-NormDepth(xImCart(3,i)))/(nPixels+1);
        end
        end
    end
    end
    end
end
imshow(depthMap2);

% Rotation = [1,0,0;0,cos(pi/12),-sin(pi/12);0,sin(pi/12),cos(pi/12)];
Rotation = [cos(pi/6) 0 sin(pi/6); 0 1 0; -sin(pi/6) 0 cos(pi/6)];
% Rotation = [cos(pi/12) -sin(pi/12) 0; sin(pi/12) cos(pi/12) 0; 0 0 1];
NewCamera = Rotation*Camera;
NewCamExtrinsic = [[CamR*inv(Rotation) CamT];0 0 0 1];
xImCartNew = projectiveCamera(CamIntrinsic,NewCamExtrinsic,P);
xImCartNew = round(xImCartNew);
xImCartNew = [xImCartNew(2,:);xImCartNew(1,:);linspace(1,size(xImCartNew,2),size(xImCartNew,2))];
xImCartNew = sortrows(xImCartNew')';
newDepth = P-repmat(NewCamera,1,(size(P,2)));
newDepth = sqrt(newDepth(1,:).^2+newDepth(2,:).^2+newDepth(3,:).^2);
minND = min(newDepth);
maxND = max(newDepth);
newNormDepth = (maxND-newDepth)./(maxND-minND);

depthMap1 = zeros(nRows,nColumns);
for i=1:size(P,2)-1
    if (xImCartNew(1,i)>=1&&xImCartNew(1,i)<=nRows&&xImCartNew(2,i)>=1&&xImCartNew(2,i)<=nColumns)
    depthMap1(xImCartNew(1,i),xImCartNew(2,i))=newNormDepth(xImCartNew(3,i));
    end
end
figure;
imshow(depthMap1);

depthMap1 = zeros(nRows,nColumns);
for i=1:size(P,2)-1
    if (xImCartNew(1,i)>=1&&xImCartNew(1,i)<=nRows&&xImCartNew(1,i+1)>=1&&xImCartNew(1,i+1)<=nRows)
    if (xImCartNew(2,i)>=1&&xImCartNew(2,i)<=nColumns&&xImCartNew(2,i+1)>=1&&xImCartNew(2,i+1)<=nColumns)
    depthMap1(xImCartNew(1,i),xImCartNew(2,i))=newNormDepth(xImCartNew(3,i));
    depthMap1(xImCartNew(1,i+1),xImCartNew(2,i+1))=newNormDepth(xImCartNew(3,i+1));
    if (xImCartNew(1,i)==xImCartNew(1,i+1))
        nPixels = xImCartNew(2,i+1)-xImCartNew(2,i)-1;
        for j=1:nPixels
            depthMap1(xImCartNew(1,i),xImCartNew(2,i)+j)=newNormDepth(xImCartNew(3,i))+...
                j*(newNormDepth(xImCartNew(3,i+1))-newNormDepth(xImCartNew(3,i)))/(nPixels+1);
        end
    end
    end
    end
end
figure;
imshow(depthMap1);


function xImCart = projectiveCamera(K,T,XCart)

%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom = [XCart; ones(1,size(XCart,2))];
%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
XRef = T*XHom;
%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
xCamHom = XRef(1:3,:)./repmat(XRef(4,:),3,1);
%TO DO move points to image coordinates xImHom by applying intrinsic matrix
xImHom = K*xCamHom;
%TO DO convert points back to Cartesian coordinates xImCart
xImCart = xImHom(1:2,:)./repmat(xImHom(3,:),2,1);
end

function mask = CreateMask(Matrix,i,j,size)
mask = zeros(size,size);
for n = 1:size 
    for m = 1:size
       mask(n,m) = Matrix(i-(size-1)/2+(n-1),j-(size-1)/2+(m-1));
    end
end
end