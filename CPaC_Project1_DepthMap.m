clear all;
close all;
 
UV_Data = load('UV_1.mat');
UV_Data = UV_Data.UV_Data;
XY_Data = load('XY_1.mat');
XY_Data = XY_Data.XY_Data;
P = load('P_1.mat');
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

depthMap2 = zeros(2048,3072);
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
    if (xImCart(1,i)>=1&&xImCart(1,i)<=2048&&xImCart(1,i+1)>=1&&xImCart(1,i+1)<=2048)
    if (xImCart(2,i)>=1&&xImCart(2,i)<=3072&&xImCart(2,i+1)>=1&&xImCart(2,i+1)<=3072)
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

depthMap1 = zeros(2048,3072);
for i=1:size(P,2)-1
    if (xImCartNew(1,i)>=1&&xImCartNew(1,i)<=2048&&xImCartNew(2,i)>=1&&xImCartNew(2,i)<=3072)
    depthMap1(xImCartNew(1,i),xImCartNew(2,i))=newNormDepth(xImCartNew(3,i));
    end
end
figure;
imshow(depthMap1);

depthMap1 = zeros(2048,3072);
for i=1:size(P,2)-1
    if (xImCartNew(1,i)>=1&&xImCartNew(1,i)<=2048&&xImCartNew(1,i+1)>=1&&xImCartNew(1,i+1)<=2048)
    if (xImCartNew(2,i)>=1&&xImCartNew(2,i)<=3072&&xImCartNew(2,i+1)>=1&&xImCartNew(2,i+1)<=3072)
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