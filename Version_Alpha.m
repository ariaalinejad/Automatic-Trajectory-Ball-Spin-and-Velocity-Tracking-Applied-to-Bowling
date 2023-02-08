clc;close all;clear all;
% Access video file
v = VideoReader('Videos/Aria_1.MOV');
%v = VideoReader('IMG_0923.MOV');
opticFlow_v = opticalFlowFarneback;
opticFlow_s = opticalFlowFarneback;


% Choose part of video
im = read(v,1);
disp('Select the region of the image you want to analyze')
[J, rect] = imcrop(im);
rect = floor(rect);
%figure;

h = figure();
movegui(h);
imshow(J);


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
initSize = 0;
% Set initial textbox values
v_ball_text = sprintf('Velocity: - [m/s]');
dir_ball_text = sprintf('Direction: - degrees');
v_spin_text = sprintf('Spin velocity: - [RPM] (- [m/s])');
dir_spin_text = sprintf('Spin Direction: - degrees');

% Parameters
MaxRadius = 60;
MinRadius = 5;
slowMotionFactor = 8;

while hasFrame(v)
    reFrame = readFrame(v);
    iter = iter+1;
    im = reFrame(rect(2):rect(2)+rect(4), rect:rect(1)+rect(3), :);
    imshow(im)
    hold on
    
    % Plot the previous spin and velocity data
    rectangle('Position',[1+initSize*2,1,size(im,2),1+initSize*2], 'FaceColor', [1 1 1],'EdgeColor','k'); 
    text(1+initSize*2, 10, v_ball_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 30, dir_ball_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 50, v_spin_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 70, dir_spin_text, 'FontSize', 8, 'Color', 'k')
    
    
    % Use hough transform to find circles (darker then the background)
    [centers, radii, metric] = imfindcircles(histeq(im),[MinRadius MaxRadius], 'ObjectPolarity','dark');
    % Keep the best circle
    [val, idx] = max(metric);
    center = centers(idx,:);
    radius = radii(idx);
    
    % Estimate the flow on the image
    flow_v = estimateFlow(opticFlow_v,rgb2gray(im));
    
    
    pause(10^-3)
    if(radius)
        % record the movement of the center and radius
        movVector(iter,1) = center(1);
        movVector(iter,2) = center(2);
        movVector(iter,3) = radius;
        
        % Find the square used for estimating the spin flow
        if (initSize == 0)
            initSize=floor(radius);
        end
        try
        imball = im(floor(center(2))-initSize:floor(center(2))+initSize, floor(center(1))-initSize:floor(center(1))+initSize, :);%im(floor(center(1))-initSize:floor(center(1))+initSize, floor(center(2))-initSize:floor(center(2))+initSize, :);
        imshow(imball); % Show region used for spin flow in upper right
        
        % Find and plot the spin flow
        flow_s = estimateFlow(opticFlow_s,rgb2gray(imball));
        plot(flow_s,'DecimationFactor',[5 5],'ScaleFactor',2);
        catch
            imball = 0;
            flow_s = 0;
        end

        %Define the rectangle around the ball
        rectForOptiFlow(1) = floor(center(1)-radius); %x init top left val
        rectForOptiFlow(2) = floor(center(2)-radius); %y init val top left
        rectForOptiFlow(3) = floor(2*radius); %how much to the right x will move
        rectForOptiFlow(4) = floor(2*radius); %How much down the y will move
        %The above code creates a rectangle around the ball with r length
        %from center at all sides. so really it is a square.
        rectangle('Position', rectForOptiFlow);
        
        
        %Make a binary mask of the ball and take out indecies
        mask = zeros([size(im,1), size(im,2)]);
        try
            for r=1:radius 
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + center(1);
                yunit = r * sin(th) + center(2);

                for i=1:length(xunit)
                    mask(round(yunit(i)),round(xunit(i))) = 1;
                end
                SE = [1 1 1; 1 1 1; 1 1 1];
                mask = imclose(mask, SE);
                [x_idx, y_idx] = find(mask == 1);
            end
        catch
            mask = ones(size(im));
        end
        
        %--------Find avg velocity vector--------------
        % Only use the flow present in the ball
        try 
        squareOrientation_v = flow_v.Orientation(x_idx', y_idx');%squareOrientation = flow.Orientation(rectForOptiFlow(2):rectForOptiFlow(2)+rectForOptiFlow(4), rectForOptiFlow(1):rectForOptiFlow(1)+rectForOptiFlow(3));
        squareMagnitude_v = flow_v.Magnitude(x_idx',y_idx');%squareMagnitude = flow.Magnitude(rectForOptiFlow(2):rectForOptiFlow(2)+rectForOptiFlow(4), rectForOptiFlow(1):rectForOptiFlow(1)+rectForOptiFlow(3));
        catch % if part of the ball is outside frame, simply use all values
            squareOrientation_v = flow_v.Orientation;
            squareMagnitude_v = flow_v.Magnitude;
        end
        
        % We find the x and y components of the average vectors 
        avg_vx = mean2(cos(squareOrientation_v).*squareMagnitude_v);
        avg_vy = mean2(sin(squareOrientation_v).*squareMagnitude_v);
        
        % Find the direction and magnitude of the average vector
        averageMagnitude_v = sqrt(avg_vx^2 + avg_vy^2);
        averageOrientation_v = atan(avg_vy/avg_vx);  
        
        % -------Calculate the velocity and direction of the ball--------
        v_ball= ((averageMagnitude_v*v.FrameRate*0.12)/radius)*slowMotionFactor;
        v_ball = v_ball / 0.2588190451; # cos(75 degrees) division
        
        % Use - to get the make the positive angular direction the same as
        % we are used to
        dir_ball = (averageOrientation_v/pi)*180;  
        
        %--------Find avg spin vector--------------
        try
        % Only use the flow present in the ball (values have to be adjusted wrt. the changed image dimentions)
        x_idx = x_idx - (floor(center(2))-initSize);
        y_idx = y_idx - (floor(center(1))-initSize);
        % In case parts of the ball are outside of the box we are looking at index is set to zero
        x_idx(x_idx < 1) = 1;x_idx(x_idx > initSize*2-1) = initSize*2+1;
        y_idx(y_idx < 1) = 1;y_idx(y_idx > initSize*2-1) = initSize*2+1;
        squareOrientation_s = flow_s.Orientation(x_idx', y_idx');
        squareMagnitude_s = flow_s.Magnitude(x_idx', y_idx');
        catch % assume error in choosing frame around ball, and set to zero
            squareOrientation_s = 0;
            squareMagnitude_s = 0;
        end
        
        % We find the x and y components of the average vectors 
        avg_sx = mean2(cos(squareOrientation_s).*squareMagnitude_s);
        avg_sy = mean2(sin(squareOrientation_s).*squareMagnitude_s);
        
        % Find the direction and magnitude of the average vector
        averageMagnitude_s = sqrt(avg_sx^2 + avg_sy^2);
        averageOrientation_s = atan(avg_sy/avg_sx);  
        
        % -------Calculate the amplitude and direction of the ball spin--------
        v_spin= ((averageMagnitude_s*v.FrameRate*0.12)/radius)*slowMotionFactor;
        rpm_spin = (v_spin*60)/(0.12*2*pi); %use these calculations later
        
        % Use - to get the make the positive angular direction the same as
        % we are used to
        dir_spin = (averageOrientation_s/pi)*180;  
        
        % ----------PLOTS-----------
        %plot of the circle
        th = 0:pi/50:2*pi;
        xunit = radius * cos(th) + center(1);
        yunit = radius * sin(th) + center(2);
        
        plot(xunit,yunit, 1,3, 'r', 'LineWidth',3);
        plot(center(1),center(2),1,3, '.', 'MarkerSize', 5);
        
        %plot velocity vector
        quiver(center(1), center(2), avg_vx*20, avg_vy*20, 'LineWidth', 2);
        %plot spin vector
        quiver(center(1), center(2), avg_sx*20, avg_sy*20, 'LineWidth', 2);
        
        % Plot a box with the magnitude and orientation of the velocity
        % vector
        v_ball_text = sprintf('Velocity: %.4f [m/s]', v_ball);
        dir_ball_text = sprintf('Direction: %.1f degrees',dir_ball);
        
        % Plot a box with the magnitude and orientation of the spin
        % vector
        v_spin_text = sprintf('Spin velocity: %.2f [RPM] (%.4f [m/s])',rpm_spin, v_spin);
        dir_spin_text = sprintf('Spin Direction: %.1f degrees',dir_spin);
    elseif(iter>1)
        % record the movement of the center and radius
        % if no circle is detected record the previous value over
        movVector(iter,1) = movVector(iter-1,1);
        movVector(iter,2) = movVector(iter-1,2);
        movVector(iter,3) = movVector(iter-1,3);
    end
        
    
     % Save the frame in structure for later saving to video file
    s(k) = getframe(h);
    k = k+1;
    %hold off
   
end
v2 = VideoReader('Videos/slowmoCut2.mp4');
figure;
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

%% HERE WE FIND PLOT THE TRAJECTORY SEEN FROM ABOVE

% Read the first video fram, uses chooses 8 points
imshow(J);

% Remove potental zeros at the start of the move vector
for i=1:size(movVector,1)
    if (movVector(i,:)==0)
        mV = movVector(i+1:end,:);
    end
end

disp('Select 8 points, then press enter');
[xi, yi] = getpts;

hold on

% plot the lines
plot([xi(1) xi(2)], [yi(1), yi(2)], 'r', 'linewidth', 5);
plot([xi(3) xi(4)], [yi(3), yi(4)], 'r', 'linewidth', 5);
plot([xi(5) xi(6)], [yi(5), yi(6)], 'b', 'linewidth', 5);
plot([xi(7) xi(8)], [yi(7), yi(8)], 'b', 'linewidth', 5);

% get the homogeneous poins
a = [xi(1); yi(1); 1];
b = [xi(2); yi(2); 1];
c = [xi(3); yi(3); 1];
d = [xi(4); yi(4); 1];
e = [xi(5); yi(5); 1];
f = [xi(6); yi(6); 1];
g = [xi(7); yi(7); 1];
h = [xi(8); yi(8); 1];

% get the two points at infinity
l1 = cross(a, b);
l2 = cross(c, d);
i = cross(l1, l2); %point at inf
l3 = cross(e, f);
l4 = cross(g, h);
j = cross(l3, l4); %point at inf
% normalize points
i = i/i(3);
j = j/j(3);

% plot lines to line at inf. and line at inf.
plot([i(1) j(1)], [i(2) j(2)], 'g--');
plot([a(1) i(1)], [a(2) i(2)], 'b');
plot([d(1) i(1)], [d(2) i(2)], 'b');
plot([f(1) j(1)], [f(2) j(2)], 'b');
plot([h(1) j(1)], [h(2) j(2)], 'b');

hold off
%%
% Find the projective transform matrix H
horizon = cross(i,j); %horizon line
horizon = horizon/horizon(3); %normalize
H = [-1 0 0; 0 -1 0; horizon(1) horizon(2) 1];
% double check that the H matrix is correct
disp('Horizon: ');
disp(vpa(inv(H)'*[0;0;1])); %vpa for higher precision
disp('Line at infinity:');
disp(vpa(H'*horizon));


% Plot the path of the ball on the original image
figure;imshow(J);
hold on

line([xi(1),xi(2)], [yi(1), yi(2)], 'Color','red','LineStyle','--', 'LineWidth', 1.5);
line([xi(3),xi(4)], [yi(3), yi(4)], 'Color','red','LineStyle','--', 'LineWidth', 1.5);
plot(mV(:,1),mV(:,2), 'Color','blue', 'LineWidth', 1.5);
legend_str{1} = 'Side line';
legend_str{2} = 'Side line';
legend_str{3} = 'Measured trajectory';
legend(legend_str);
hold off

% find the the new path and lines seen from above
a = H*[xi(1); yi(1); 1];
b = H*[xi(2); yi(2); 1];
c = H*[xi(3); yi(3); 1];
d = H*[xi(4); yi(4); 1];
a = a/a(3); b = b/b(3); c = c/c(3); d = d/d(3); % normalize
movVector_new = mV; 
movVector_new(:,3) = 1;
movVector_new = (H*(movVector_new'))';
movVector_new = movVector_new./movVector_new(:,3);

% Plot the lines and path seen from above
figure;
hold on
line([a(1), b(1)], [a(2), b(2)], 'Color', 'blue');
line([c(1), d(1)], [c(2), d(2)], 'Color', 'blue');
plot(movVector_new(:,1),movVector_new(:,2), 'Color', 'red');
legend_str{1} = 'Side line';
legend_str{2} = 'Side line';
legend_str{3} = 'Measured trajectory';
legend(legend_str);
hold off
