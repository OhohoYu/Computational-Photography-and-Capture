close all;
clear all;

% im1 = imread('./own1.JPG'); 
% im1 = imresize(im1,0.5);
% im1 = im2double(im1);
% 
% im2 = imread('./own2.JPG'); 
% im2 = imresize(im2,0.5);
% im2 = im2double(im2);
% 
% im3 = imread('./own3.JPG'); 
% im3 = imresize(im3,0.5);
% im3 = im2double(im3);
% 
% H_p1_c = Homography_Pi_c(im1);
% H_p2_c = Homography_Pi_c(im2);
% H_p3_c = Homography_Pi_c(im3);

%*******************************************************************
% Test Homographies
% Im_P_p_projected = zeros(768,1024,3);
% 
% for cY=1:imY
%     for cX=1:imX
%         HomP1 = [cX;cY;1];
%         HomP2 = H_p1_c*HomP1;
%         CartP2 = HomP2(1:2,:)/HomP2(3,:);
%         CartP2 = round(CartP2);
%         CartP2(CartP2<1) = 1;
%         Im_P_p_projected(CartP2(2),CartP2(1),:) = im1(cY,cX,:);
%     end
% end
% 
% imshow(Im_P_p_projected);
%*******************************************************************   

H_p1_c = load('own_H_p1_c.mat');
H_p1_c = H_p1_c.H_p1_c;
H_p2_c = load('own_H_p2_c.mat');
H_p2_c = H_p2_c.H_p2_c;
H_p3_c = load('own_H_p3_c.mat');
H_p3_c = H_p3_c.H_p3_c;

H_p2_p1 = H_p2_c*inv(H_p1_c);
H_p3_p1 = H_p3_c*inv(H_p1_c);


P_c_printed1 = imread('./own1.JPG'); 
P_c_printed1 = imresize(P_c_printed1,0.5);
P_c_printed1 = im2double(P_c_printed1);
P_p1_printed = CreateImage(P_c_printed1, H_p1_c);
imshow(P_p1_printed);
imwrite(P_p1_printed,'ownVirtualPattern1.jpg');

P_c_printed2 = imread('./own2.JPG'); 
P_c_printed2 = imresize(P_c_printed2,0.5);
P_c_printed2 = im2double(P_c_printed2);
P_p2_printed = CreateImage(P_c_printed2, H_p2_c);
figure;
imshow(P_p2_printed);
imwrite(P_p2_printed,'ownVirtualPattern2.jpg');

P_c_printed3 = imread('./own3.JPG'); 
P_c_printed3 = imresize(P_c_printed3,0.5);
P_c_printed3 = im2double(P_c_printed3);
P_p3_printed = CreateImage(P_c_printed3, H_p3_c);
figure;
imshow(P_p3_printed);
imwrite(P_p3_printed,'ownVirtualPattern3.jpg');


%==========================================================================
function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding matchin in 
%pts2Cart

%**** TO DO ****;
%first turn points to homogeneous
%then construct A matrix which should be (10 x 9) in size
nData = size(pts1Cart,2);
for i=1:nData
    A(2*i-1,:) = [pts1Cart(1,i) pts1Cart(2,i) 1 0 0 0 -pts1Cart(1,i)*pts2Cart(1,i) -pts1Cart(2,i)*pts2Cart(1,i) -pts2Cart(1,i)];
    A(2*i,:) = [0 0 0 pts1Cart(1,i) pts1Cart(2,i) 1 -pts1Cart(1,i)*pts2Cart(2,i) -pts1Cart(2,i)*pts2Cart(2,i) -pts2Cart(2,i)];    
    
end
%solve Ah = 0 by calling
%h = solveAXEqualsZero(A); (you have to write this routine too - see below)
h = solveAXEqualsZero(A);
%reshape h into the matrix H
H = reshape(h,[3 3])';
%Beware - when you reshape the (9x1) vector x to the (3x3) shape of a homography, you must make
%sure that it is reshaped with the values going first into the rows.  This
%is not the way that the matlab command reshape works - it goes columns
%first.  In order to resolve this, you can reshape and then take the
%transpose
end

%==========================================================================
function x = solveAXEqualsZero(A);
[U,L,V] = svd(A);
x = V(:,end);
end

%==========================================================================

function H_p1_c = Homography_Pi_c(im1)
imshow(im1);
[x1 y1] = ginput(4);
% P-c-projected
P1 = [x1';y1'];

% P-p-projected: 1280*720
% Order: Top-Left, Top-Right, Bot-Right, Bot-Left
P2 = [69.3333,437.6666,437.6666,69.3333;187.3333,187.3333,528.6666,528.6666];

H_p1_c = calcBestHomography(P1,P2);
%H_c_p1 = calcBestHomography(P2,P1);
end

%==========================================================================
function P_p_printed = CreateImage(P_c_printed, H_p_c)
[imY, imX, nDim] = size(P_c_printed);
P_p_printed = zeros(720,1280,3);
for cY=1:imY
    for cX=1:imX
        HomP1 = [cX;cY;1];
        HomP2 = H_p_c*HomP1;
        CartP2 = HomP2(1:2,:)/HomP2(3,:);
        CartP2 = round(CartP2);
        if CartP2(1)>1280
            CartP2(1) = 1280;
        end
        if CartP2(2)>720
            CartP2(2) = 720;
        end
        CartP2(CartP2<1) = 1;
        P_p_printed(CartP2(2),CartP2(1),:) = P_c_printed(cY,cX,:);
    end
end
end