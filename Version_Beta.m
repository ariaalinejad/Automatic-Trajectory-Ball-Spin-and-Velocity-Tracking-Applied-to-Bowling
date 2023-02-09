clc;close all;clear all;
% Access video file
%v = VideoReader('Videos/IMG_0923.MOV');
v = VideoReader('Videos/Aria_1.MOV');

% Choose frame of video
im = read(v,v.NumFrames/2);
imshow(im)
title('Select the region of the image you want to analyze, then double kick on it');
hold on
[J, rect] = imcrop(im);
rect = floor(rect);
hold off
%------Load values during testing to save time----
%rect = load('R/vars.mat', 'rect').rect;
%J = im(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3), :);
%-----------------------------------------------------------------

imshow(J);

%---------Select region we are interested in (the lane)----------
title('Select the four angles of the lane, then press enter');
[xi, yi] = getpts;
%-----Load values during testing to save time----
%yi = load('R/vars.mat', 'yi').yi;
%xi = load('R/vars.mat', 'xi').xi;
%-----------------------------------------------------------------
hold on

lane_x = [xi(1), xi(2), xi(3), xi(4)];
lane_y = [yi(1), yi(2), yi(3), yi(4)];

U = roipoly(J,lane_x,lane_y);
U = bwareafilt(U,1); 
U = bwconvhull(U); 
croppedImg = immultiply(J,repmat(U,[1 1 3])); 
imshow(croppedImg);

%%

% Define figure and starting frame
h = figure();
movegui(h);
im = read(v,1);

% Preallocate structure to store video frames
% Note that this is just an initial guess for preallocation based on
% duration and framerate, the video may have fewer or more frames
nFrames = ceil(v.FrameRate*v.Duration);
s(nFrames) = struct('cdata',[],'colormap',[]);


% Loop through video, grabbing frames and updating plots
k = 1;
iter = 0;
movVector = zeros(nFrames,3); % Vector storing centers and radii
rectBall = [0,0,0,0]; % rectangle plotted around ball
initSize = 0; % radius of the first detected all
% Set initial textbox values
v_ball_text = sprintf('Velocity: - [m/s]');
dir_ball_text = sprintf('Direction: - degrees');
v_spin_text = sprintf('Spin velocity: - [RPM] (- [m/s])');
dir_spin_text = sprintf('Spin Direction: - degrees');
% Define Optical flow variables
opticFlow_v = opticalFlowFarneback; % Optical flow for velocity
opticFlow_s = opticalFlowFarneback; % Optical flow for spin

% Parameters
MaxRadius = 60; % Found from testing, radius of biggest ball
MinRadius = 7; % Found from testing, radius of smallest ball
slowMotionFactor = 8; % Found from knowledge og original video frame rate 

while hasFrame(v)
    reFrame = readFrame(v);
    iter = iter+1;
    im = reFrame(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3), :);
    imshow(im)
    hold on
    
    % Only look at area of image where the lane is
    imCrop = immultiply(im,repmat(U,[1 1 3]));
    
    % Plot the previous spin and velocity data (so that the values will flicker less)
    rectangle('Position',[1+initSize*2,1,size(im,2),1+initSize*2], 'FaceColor', [1 1 1],'EdgeColor','k'); 
    text(1+initSize*2, 10, v_ball_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 30, dir_ball_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 50, v_spin_text, 'FontSize', 8, 'Color', 'k')
    text(1+initSize*2, 70, dir_spin_text, 'FontSize', 8, 'Color', 'k')
    
    % Use hough transform to find circles (darker then the background)
    % We use adapthisteq on the image for better results
    if (iter>1)
        [centers, radii, metric] = imfindcircles(adapthisteq(rgb2gray(imCrop)),[MinRadius MaxRadius], 'ObjectPolarity','dark');
        % Keep circle most similar to last circle
        [val, idx] = min(abs(radii-movVector(iter-1,3)));
        center = centers(idx,:);
        radius = radii(idx);
    else
        [centers, radii, metric] = imfindcircles(adapthisteq(rgb2gray(imCrop)),[MaxRadius*0.8 MaxRadius], 'ObjectPolarity','dark');
        % Keep the best circle
        [val, idx] = max(metric);
        center = centers(idx,:);
        radius = radii(idx);
    end
   
    
    % Estimate the flow on the image
    flow_v = estimateFlow(opticFlow_v,rgb2gray(imCrop));
    
    pause(10^-3)
    if(radius) % checks if we found a circle or not
        % record the movement of the center and radius
        movVector(iter,1) = center(1);
        movVector(iter,2) = center(2);
        movVector(iter,3) = radius;
        
        % Find the square used for estimating the spin flow
        if (initSize == 0)
            initSize=floor(radius);
        end
        
        % use try catch to prevent error when part of the ball is out of frame
        % imball is the image around the ball used for spin detection
        try 
            imball = im(floor(center(2))-initSize:floor(center(2))+initSize, floor(center(1))-initSize:floor(center(1))+initSize, :);%im(floor(center(1))-initSize:floor(center(1))+initSize, floor(center(2))-initSize:floor(center(2))+initSize, :);
            imshow(imball); % Show region used for spin flow in upper right

            % Find and plot the spin flow
            flow_s = estimateFlow(opticFlow_s,rgb2gray(imball));
            plot(flow_s,'DecimationFactor',[5 5],'ScaleFactor',2);
        catch
            flow_s = 0;
        end

        %Define the rectangle around the ball
        rectBall(1) = floor(center(1)-radius); %x init top left val
        rectBall(2) = floor(center(2)-radius); %y init val top left
        rectBall(3) = floor(2*radius); %how much to the right x will move
        rectBall(4) = floor(2*radius); %How much down the y will move
        %The above code creates a rectangle around the ball with r length
        %from center at all sides. so really it is a square.
        rectangle('Position', rectBall); %plots the rectangle
        
        
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
        % try catch to prevent error when part of ball if out of frame
        try 
            Vx_v = flow_v.Vx(x_idx', y_idx');
            Vy_v = flow_v.Vy(x_idx', y_idx');
        catch % if part of the ball is outside frame, simply use all values
            Vx_v = flow_v.Vx;
            Vy_v = flow_v.Vy;
        end
        
        
        % We find the x and y components of the average vectors 
        avg_vx = mean2(Vx_v);
        avg_vy = mean2(Vy_v);
        
        % Find the direction and magnitude of the average vector
        averageMagnitude_v = sqrt(avg_vx^2 + avg_vy^2);
        
        % Change axis so that they are as we are used to
        if(avg_vx>0 && avg_vy>0)
            averageOrientation_v =  atan(avg_vy/avg_vx);
        elseif (avg_vx<0 && avg_vy>0)
            averageOrientation_v =  pi - atan(avg_vy/(-avg_vx));
        elseif (avg_vx<0 && avg_vy<0)
            averageOrientation_v =   pi + atan((-avg_vy)/(-avg_vx));
        elseif (avg_vx>0 && avg_vy<0)
            averageOrientation_v =  2*pi - atan((-avg_vy)/(avg_vx));
        end
        averageOrientation_v = 2*pi - averageOrientation_v;
        
        
        % -------Calculate the velocity and direction of the ball--------
        v_ball= ((averageMagnitude_v*v.FrameRate*0.12)/radius)*slowMotionFactor;
        v_ball = v_ball/0.2588190451; % cos(75 degrees) division
        dir_ball = (averageOrientation_v/pi)*180; 
        
        %--------Find avg spin vector--------------
        try
            % Only use the flow present in the ball (values have to be adjusted wrt. the changed image dimentions)
            x_idx = x_idx - (floor(center(2))-initSize);
            y_idx = y_idx - (floor(center(1))-initSize);
            % In case parts of the ball are outside of the box we are
            % looking at index is set to one or end value
            x_idx(x_idx < 1) = 1;x_idx(x_idx > initSize*2-1) = initSize*2+1;
            y_idx(y_idx < 1) = 1;y_idx(y_idx > initSize*2-1) = initSize*2+1;
            % Find flow only at binary ball mask
            Vx_s= flow_s.Vx(x_idx', y_idx');
            Vy_s = flow_s.Vy(x_idx', y_idx');
            M_s = flow_s.Magnitude(x_idx', y_idx');
        catch % assume error in choosing frame around ball, and set to zero
            Vx_s = 0;
            Vy_s = 0;
            M_s = 0;
        end
        
        % Only use the highest (most prominant spin values in the calculation)
        M_ss = maxk(M_s, 10, 1);
        [~, m_I] = maxk(M_ss, 10, 2);
        Vx_s_m = Vx_s(m_I);
        Vy_s_m = Vy_s(m_I);
        
        % We find the x and y components of the average vectors 
        avg_sx = mean2(Vx_s_m);
        avg_sy = mean2(Vy_s_m);
        
        % Find the direction and magnitude of the average vector
        averageMagnitude_s = sqrt(avg_sx^2 + avg_sy^2);
        
        % Change axis so that they are as we are used to
        if(avg_sx>0 && avg_sy>0)
            averageOrientation_s =  atan(avg_sy/avg_sx);
        elseif (avg_sx<0 && avg_sy>0)
            averageOrientation_s =  pi - atan(avg_sy/(-avg_sx));
        elseif (avg_sx<0 && avg_sy<0)
            averageOrientation_s =   pi + atan((-avg_sy)/(-avg_sx));
        elseif (avg_sx>0 && avg_sy<0)
            averageOrientation_s =  2*pi - atan((-avg_sy)/(avg_sx));
        else
            averageOrientation_s = 0;
        end
        averageOrientation_s = 2*pi - averageOrientation_s;
         
        
        % -------Calculate the amplitude and direction of the ball spin--------
        v_spin= ((averageMagnitude_s*v.FrameRate*0.12)/radius)*slowMotionFactor;
        rpm_spin = (v_spin*60)/(0.12*2*pi); 
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
    hold off
   
end

v2 = VideoReader('Videos/Aria_1.MOV');
v2 = readFrame(v2);
figure;
imshow(v2(rect(2):rect(2)+rect(4), rect:rect(1)+rect(3), :))
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


%% HERE WE FIND AND PLOT THE TRAJECTORY SEEN FROM ABOVE

% Read a single video frame
imshow(J);

%-------To allow for testing without running the whole code------
%movVector = load('R/vars.mat', 'movVector').movVector;
%----------------------------------------------------------------

% Remove potental zeros
mV = reshape(movVector(movVector>0), [size(movVector(movVector>0),1)/3, 3]);

%title('Select points (eg. 32 pointsthat form 8 right angles), first 4 must be around lane edges!, then press enter');
%[xii, yii] = getpts;
%---To save you the hastle of choosing all of the points, here you can load
% them (given that the same video crop is used)
xii = load('R/vars.mat', 'xii').xii;
yii = load('R/vars.mat', 'yii').yii;
%-------------------------------------------------------------------
hold on

% Initialize
N = size(xii,1); % num points
a = zeros(N,3); % points in homogeneous coordinates
l = zeros(3,N/4); % lines l
m = zeros(3,N/4); % lines m (orthogonal to l)
A = zeros(N/4,6); 

for i=1:N
    a(i,:) = [xii(i); yii(i); 1];
end
for i=1:2:(N/2-1)
    plot([xii(i) xii(i+1)], [yii(i), yii(i+1)], 'linewidth', 5);
end
for i=1:4:N
    l(:,i) = cross(a(i,:)', a(i+1,:)');
    m(:,i) = cross(a(i+2,:)', a(i+3,:)');
end
for i=1:N/4
    A(i, :) = [l(1,i)*m(1,i),0.5*(l(1,i)*m(2,i)+l(2,i)*m(1,i)),l(2,i)*m(2,i),...
            0.5*(l(1,i)*m(3,i)+l(3,i)*m(1,i)),  0.5*(l(2,i)*m(3,i)+l(3,i)*m(2,i)), l(3,i)*m(3,i)];
end
         
%%

[~,~,v] = svd(A); 
sol = v(:,end); 
imDCCP = [sol(1)  , sol(2)/2, sol(4)/2;...
    sol(2)/2, sol(3)  , sol(5)/2;...
    sol(4)/2, sol(5)/2  sol(6)];


[U,D,V] = svd(imDCCP);
D(3,3) = 1;
A = U*sqrt(D);

C = [eye(2),zeros(2,1);zeros(1,3)];
min(norm(A*C*A' - imDCCP),norm(A*C*A' + imDCCP))

H = inv(A);
min(norm(H*imDCCP*H'./norm(H*imDCCP*H') - C./norm(C)),norm(H*imDCCP*H'./norm(H*imDCCP*H') + C./norm(C)))


tform = projective2d(H');
K = imwarp(J,tform);

figure;
imshow(K);

hold off
%%

b = zeros(3,8); % transformed points

for i=1:size(b,2)
    b(:,i) = H*a(i,:)';
    b(:,i) = b(:,i)/b(3,i); % normalize
end

figure;imshow(J);
hold on
for i=1:2:size(b,2) % assume 4 first selected lines as lines around lane
    line([a(i,1),a(i+1,1)], [a(i,2),a(i+1,2)], 'LineWidth', 3, 'Color', 'blue');
end
plot(mV(:,1),mV(:,2), 'Color','red', 'LineWidth', 1.5); % plot trajectory on lane image
hold off

% Plot the lines and path of original video
figure;
subplot(1,2,1);
hold on 
for i=1:2:size(b,2)
    line([a(i,1),a(i+1,1)], [a(i,2),a(i+1,2)], 'Color', 'blue');
end
plot(mV(:,1),mV(:,2), 'Color','red');
set(gca, 'YDir','reverse'); % We invert the axis since a normal image has y axis inverted when shown
hold off

% Find transformed trajectory
mV_new = mV; 
mV_new(:,3) = 1;
mV_new = (H*mV_new')';
mV_new = mV_new./mV_new(:,3);

% Plot transfromed lines and trajectory
b = -b; % The values are set to negative just to make the plot more plesant to campare
mV_new = -mV_new;
subplot(1,2,2);
hold on
for i=1:2:size(b,2)
    line([b(1,i),b(1,i+1)], [b(2,i),b(2,i+1)], 'Color', 'blue');
end
plot(mV_new(:,1),mV_new(:,2), 'Color', 'red');
set(gca, 'YDir','reverse');
hold off