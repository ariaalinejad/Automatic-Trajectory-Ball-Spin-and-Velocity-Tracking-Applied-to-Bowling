%%Clear all
close all;clear;clc;

%% Reading the video
video = vision.VideoFileReader('Videos/slowmoCut.mp4');
video_frames = read(VideoReader('Videos/slowmoCut.mp4'));
%Now video_frames contains the frames we need.

player = vision.DeployableVideoPlayer('Location', [10,100]); %Video player


%%Detect the ball
fgDetector = vision.ForegroundDetector('NumTrainingFrames', 10, 'InitialVariance', 0.0015); %0.05
fgPlayer = vision.DeployableVideoPlayer('Location', player.Location + [500 120]);
blobAnalyzer = vision.BlobAnalysis('AreaOutputPort', true, 'MinimumBlobArea', 10, 'CentroidOutputPort', true);


frameNr = 1; %starts the count on framenumber on 1 
crop_x = 650;
crop_y = 350;
crop_width = 200;
crop_height = 300;
counter = 0;
while ~isDone(video)
    
    im = step(video);
    image = imcrop(im,[crop_x, crop_y, crop_width, crop_height]);
    I = rgb2gray(image);
    
    
    %detect ball
    str = strel("disk",6);
    fgMask = step(fgDetector,I); %masking image to get at binary image
    fgMask = bwareaopen(fgMask,100); %removes small objects from binary image, e.g. noise in background
    fgMask = imopen(fgMask,str);
    [~, detection] = step(blobAnalyzer,fgMask); %detection is the values of detected objected
    step(fgPlayer,fgMask);

    if ~isempty(detection) %if there is a detected object in the dataframe
        counter = counter + 1;
        position = detection(1,:);
        pos(frameNr,:) = detection(1,:);
        position(:,3) = 10;
        combinedImage = insertShape(image,'circle',position); %,'Ball'
        step(player,combinedImage);
    else
        step(player, image);
        pos(frameNr,:) = [NaN NaN];
    end
    step(fgPlayer,fgMask);
    
    step(player,image);
    
    frameNr = frameNr + 1; %increase the index/framenumber after each iteration
end
%% Visualize the trajectory
figure, imshow(image)
hold on
plot(pos(:,1),pos(:,2),'r-o','LineWidth',2)
legend('Measured trajectory')

        