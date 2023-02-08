clc;close all;clear all; 

v = VideoReader('Videos/Aria_1.MOV');
%v = VideoReader('IMG_0923.MOV');
im = read(v,500);
J = im(170:170+860, 600:600+709, :);
figure;imshow(J);
hold all

disp('Select 8 points, then press enter');
[xi, yi] = getpts;

%%
v = VideoReader('Videos/Aria_1.MOV');
im = read(v,500);
J = im(170:170+860, 600:600+709, :);
figure;imshow(J);
hold all

N = size(xi,1);
a = zeros(N,3);
l = zeros(3,N/4);
m = zeros(3,N/4);
A = zeros(N/4,6);

for i=1:N
    a(i,:) = [xi(i); yi(i); 1];
end
for i=1:2:(N/2-1)
    plot([xi(i) xi(i+1)], [yi(i), yi(i+1)], 'linewidth', 5);
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
            
[~,~,v] = svd(A); %
sol = v(:,end); %sol = (a,b,c,d,e,f)  [a,b/2,d/2; b/2,c,e/2; d/2 e/2 f];
imDCCP = [sol(1)  , sol(2)/2, sol(4)/2;...
    sol(2)/2, sol(3)  , sol(5)/2;...
    sol(4)/2, sol(5)/2  sol(6)];


[U,D,V] = svd(imDCCP);
D(3,3) = 1;
A = U*sqrt(D);

C = [eye(2),zeros(2,1);zeros(1,3)];
min(norm(A*C*A' - imDCCP),norm(A*C*A' + imDCCP))

H = inv(A); % rectifying homography
min(norm(H*imDCCP*H'./norm(H*imDCCP*H') - C./norm(C)),norm(H*imDCCP*H'./norm(H*imDCCP*H') + C./norm(C)))


tform = projective2d(H');
K = imwarp(J,tform);

figure;
imshow(K);

hold off
%%
b = zeros(3,8);

for i=1:size(b,2)
    b(:,i) = H*a(i,:)';
    b(:,i) = b(:,i)/b(3,i);
end

figure;imshow(J);
hold on
for i=1:2:size(b,2)
    line([a(i,1),a(i+1,1)], [a(i,2),a(i+1,2)], 'Color', 'blue');
end
hold off

figure;
subplot(1,2,1);
for i=1:2:size(b,2)
    line([b(1,i),b(1,i+1)], [b(2,i),b(2,i+1)], 'Color', 'blue');
end
hold on 
subplot(1,2,2);
for i=1:2:size(b,2)
    line([a(i,1),a(i+1,1)], [a(i,2),a(i+1,2)], 'Color', 'blue');
end
hold off
