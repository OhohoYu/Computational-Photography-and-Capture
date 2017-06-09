clear all;
close all;
UV_Data = load('own_UV.mat');
UV_Data = UV_Data.UV_Data;
XY_Data = load('own_XY.mat');
XY_Data = XY_Data.XY_Data;
P = load('own_P.mat');
P = P.P;

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

nRows = 1728;
nColumns = 2592;
depthMap2 = zeros(nRows,nColumns);

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