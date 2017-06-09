clear all;
close all;
fileName = 'DroneMov.mp4';
obj = VideoReader(fileName);
numFrames = obj.NumberOfFrames;
FrameRate = obj.FrameRate;
% CompressedTime = 15;
% SamplingInterval = floor(obj.Duration/CompressedTime);
% i=0;
% for k = 1 :14:numFrames
%     fprintf('processing No.%d \n',k);
%     frame = read(obj,k);
%     frame=imresize(frame,0.3);
%     %imshow(frame);
%     if i<10
%     imwrite(frame,strcat('testFootage00',num2str(i),'.jpg'),'jpg');
%     end
%     if i<100&&i>=10
%     imwrite(frame,strcat('testFootage0',num2str(i),'.jpg'),'jpg');
%     end
%     if i<1000&&i>=100
%     imwrite(frame,strcat('testFootage',num2str(i),'.jpg'),'jpg');
%     end
% %     if i<10000&&i>=1000
% %     imwrite(frame,strcat('testFootage',num2str(i),'.jpg'),'jpg');
% %     end
% %     if i<100000&&i>=10000
% %     imwrite(frame,strcat('testFootage',num2str(i),'.jpg'),'jpg');
% %     end
%     i=i+1;
% end

movOriginal = load_sequence_color('./Drones','testFootage', 0, 833, 3, 'jpg');
[nRows, nColumns, nColors, nFrames] = size(movOriginal);


movGray = zeros(nRows,nColumns,nFrames);
for i=1:nFrames
    movGray(:,:,i) = rgb2gray(movOriginal(:,:,:,i));
end