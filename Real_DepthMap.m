UV_Data = load('crayon4_UV.mat');
UV_Data = UV_Data.UV_Data;
XY_Data = load('crayon4_XY.mat');
XY_Data = XY_Data.XY_Data;
P = load('crayon4_P.mat');
P = P.P;

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

nRows = 3456;
nColumns = 5184;
depthMap2 = zeros(nRows,nColumns);

xImCart = [xImCart(2,:);xImCart(1,:);linspace(1,size(xImCart,2),size(xImCart,2))];
xImCart = sortrows(xImCart')';
for i=1:size(P,2)-1
    if (xImCart(1,i)>=1&&xImCart(1,i)<=nRows&&xImCart(1,i+1)>=1&&xImCart(1,i+1)<=nRows)
    if (xImCart(2,i)>=1&&xImCart(2,i)<=nColumns&&xImCart(2,i+1)>=1&&xImCart(2,i+1)<=nColumns)
    
    depthMap2(xImCart(1,i),xImCart(2,i))=NormDepth(xImCart(3,i));
    depthMap2(xImCart(1,i+1),xImCart(2,i+1))=NormDepth(xImCart(3,i+1));
    if (xImCart(1,i)==xImCart(1,i+1))
        if (abs(xImCart(2,i+1)-xImCart(2,i))<=100)
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