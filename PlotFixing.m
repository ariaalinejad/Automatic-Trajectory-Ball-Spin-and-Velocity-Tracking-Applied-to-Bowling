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
optiFlowVec = zeros(nFrames,2); % One vector is orientation, the other is magnitude
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

    frameGray = rgb2gray(frameRGB); 
    
    flow = estimateFlow(opticFlow,frameGray);
    
    imshow(im)
    hold on
    
    %Plot flow of whole video frame
    %plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2);
    
    pause(10^-3)
    %viscircles(center, radius,'EdgeColor','b');
    if(radius)
        
        th = 0:pi/50:2*pi;
        xunit = radius * cos(th) + center(1);
        yunit = radius * sin(th) + center(2);

        %Define the rectangle in which we want to find the avg optical flow

        rectForOptiFlow(1) = floor(center(1)-radius); %x init top left val
        rectForOptiFlow(2) = floor(center(2)-radius); %y init val top left
        rectForOptiFlow(3) = floor(2*radius); %how much to the right x will move
        rectForOptiFlow(4) = floor(2*radius); %How much down the y will move
        
        %The above code creates a rectangle around the ball with r length
        %from center at all sides. so really it is a square.

        rectangle('Position', rectForOptiFlow);

        squareOrientation = flow.Orientation(rectForOptiFlow(2):rectForOptiFlow(2)+rectForOptiFlow(4), rectForOptiFlow(1):rectForOptiFlow(1)+rectForOptiFlow(3));

        % Calculate the average
        % orientations goes from -pi to pi
        averageOrientation = mean2(squareOrientation);

        squareMagnitude = flow.Magnitude(rectForOptiFlow(2):rectForOptiFlow(2)+rectForOptiFlow(4), rectForOptiFlow(1):rectForOptiFlow(1)+rectForOptiFlow(3));
        
        % Calculate the average
        % calculates so that vectors in opposite directons cancel
        averageMagnitude = mean2(cos(squareOrientation).*squareMagnitude);
        
        % We find the x and y magnitudes of the spin vector. We mutiply
        % with 50 to get a more visible vector on the image. We use minus
        % as these calculations are based on a standard coordinate system,
        % but the coordinate system of the image has zero in top left and
        % x/y increasing as we go down and to the right
        vx = -averageMagnitude*cos(averageOrientation)*50;
        vy = -averageMagnitude*cos(averageOrientation)*50;
        
        %plot spin vector
        quiver(center(1), center(2), vx, vy, 'LineWidth', 2);
        
        %plot of the circle
        plot(xunit,yunit, 1,3, 'r', 'LineWidth',3);
        plot(center(1),center(2),1,3, '.', 'MarkerSize', 5);
        

        
    end
    
    
    
    
    
    
    
     % Save the frame in structure for later saving to video file
    s(k) = getframe(h);
    k = k+1;
    hold off
   
end
v2 = VideoReader('Videos/slowmoCut2.mp4');
figure
imshow(readFrame(v2))
hold on
plot(movVector(:,1),movVector(:,2),'r-o','LineWidth',2);
legend('Measured trajectory');
hold off
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



% Need to estimate the spin from the flow variable. 
% Idea: Look in the area inside the circle radius
% Use trigonometry of some sort to find the spin for a certain point in the
% circle. Maybe we just need to sum all the magnitude and orientation
% vectors.

% De andre tar bare average orientation og magnitude i sirkelområdet. så da
% må vi bare gjøre det.
% Bør gjøre per frame.
