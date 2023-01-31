clc;close all;clear all;
% Access video file
v = VideoReader('Videos/slowmoCut2.mp4');
opticFlow = opticalFlowFarneback;
h = figure();
movegui(h);




% Preallocate structure to store video frames
% Note that this is just an initial guess for preallocation based on
% duration and framerate, the video may have fewer or more frames
nFrames = ceil(v.FrameRate*v.Duration);
s(nFrames) = struct('cdata',[],'colormap',[]);


% Loop through video, grabbing frames and updating plots
k = 1;
movVector = zeros(nFrames,3);
iter = 0;
rectForOptiFlow = [0,0,0,0];
while hasFrame(v)
    reFrame = readFrame(v);
    iter = iter+1;
    im = reFrame;
    im = histeq(im);

    [centers, radii, metric] = imfindcircles(im,[10 150], 'ObjectPolarity','dark');%imfindcircles(I,[15 150], 'ObjectPolarity','dark');

    [val, idx] = max(metric);
    center = centers(idx,:);
    radius = radii(idx);
    movVector(iter,1) = center(1);
    movVector(iter,2) = center(2);
    movVector(iter,3) = radius;

    frameRGB = reFrame;

    frameGray = im2gray(frameRGB); 
    
    flow = estimateFlow(opticFlow,frameGray);
    imshow(im)
    hold on
    
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2);
    
    pause(10^-3)
    %viscircles(center, radius,'EdgeColor','b');
    if(radius)
        
        th = 0:pi/50:2*pi;
        xunit = radius * cos(th) + center(1);
        yunit = radius * sin(th) + center(2);
        plot(xunit,yunit, 1,3, 'r', 'LineWidth',3);
        plot(center(1),center(2),1,3, '.', 'MarkerSize', 5);
    end
     % Save the frame in structure for later saving to video file
    s(k) = getframe(h);
    k = k+1;
    hold off
   

    

   

end

plot(movVector(:,1),movVector(:,2),'r-o','LineWidth',2);
legend('Measured trajectory');
% Remove any unused structure array elements
s(k:end) = [];

% Open a new figure and play the movie from the structure
hFig2 = figure;
movie(hFig2,s,1,v.FrameRate);

% Write to a video file
% This could be done within the original loop, but I wanted to show it
% separately
vOut = VideoWriter('Videos/slowmoCut3','MPEG-4');
vOut.FrameRate = v.FrameRate;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)
