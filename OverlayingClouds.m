close all;
clear all;
UV_Data = load('notebook_UV.mat');
UV_Data = UV_Data.UV_Data;
XY_Data = load('notebook_XY.mat');
XY_Data = XY_Data.XY_Data;
P = load('notebook_P.mat');
P = P.P;

CamIntrinsic = [4800.90780962525,0,1535.60950558134;0,4800.82059379274,1022.01737766108;0,0,1];
ProjIntrinsic = [3366.73320708613,0,511.500000000000;0,3246.86206530026,383.500000000000;0,0,1];
CamT = [ -154.928480 ;	 -155.396211 ;	 1103.799824 ]/3000;
ProjT = [ -72.123869 ;	 -35.869205 ;	 1111.176760 ]/1000;
CamR = [ -0.002008 	 0.999981 	 0.005888;
         0.946347 	 0.003803 	 -0.323130;
        -0.323146 	 0.004923 	 -0.946336 ];
ProjR = [ 0.429563 	 0.719805 	 0.545304;
            0.902282 	 -0.366800 	 -0.226595;
        0.036913 	 0.589355 	 -0.807030 ];

% CamT = [-0.10157115; -0.10139455; 0.49625999];
% ProjT = [-0.14915472; -0.06104425;  1.36014771];
% CamR = [0.99998617,      -0.00475739,       0.00223672;
%         0.00382316,       0.95016861,       0.31171292;
%        -0.00360820,      -0.31170005,       0.95017368];
% ProjR = [0.72119248,       0.44233182,      -0.53312665;
%         -0.36164442,       0.89680630,       0.25485638;
%          0.59084243,       0.00900178,       0.80673677];
    
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
P2 = w;
scatter3(P(1,:),P(2,:),P(3,:));
hold on;
scatter3(P2(1,:),P2(2,:),P2(3,:),'r');

Camera = -CamR'*CamT;
CamExtrinsic = [[CamR CamT];0 0 0 1];
ProjExtrinsic = [[ProjR ProjT];0 0 0 1];
xImCart = projectiveCamera(CamIntrinsic,CamExtrinsic,P2);
xImCart = round(xImCart);
Depth = P2-repmat(Camera,1,(size(P2,2)));
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
for i=1:size(P2,2)-1
    if (xImCart(1,i)>=1&&xImCart(1,i)<=nRows&&xImCart(1,i+1)>=1&&xImCart(1,i+1)<=nRows)
    if (xImCart(2,i)>=1&&xImCart(2,i)<=nColumns&&xImCart(2,i+1)>=1&&xImCart(2,i+1)<=nColumns)
    depthMap2(xImCart(1,i),xImCart(2,i))=NormDepth(xImCart(3,i));
    depthMap2(xImCart(1,i+1),xImCart(2,i+1))=NormDepth(xImCart(3,i+1));
    if (xImCart(1,i)==xImCart(1,i+1))
        nPixels = xImCart(2,i+1)-xImCart(2,i)-1;
        for j=1:nPixels
            depthMap2(xImCart(1,i),xImCart(2,i)+j)=NormDepth(xImCart(3,i))+...
                j*(NormDepth(xImCart(3,i+1))-NormDepth(xImCart(3,i)))/(nPixels+1);
        end
    end
    end
    end
end
figure;
imshow(depthMap2);

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