%% setup
close all; clear all; clc;

% Load the camera data
load('cam1_1.mat');
load('cam2_1.mat');
load('cam3_1.mat');

[m1,n1]=size(vidFrames1_1(:,:,1,1)); % compute data size
[m2,n2]=size(vidFrames2_1(:,:,1,1));
[m3,n3]=size(vidFrames3_1(:,:,1,1));

numf1 = size(vidFrames1_1,4); % cam1
numf2 = size(vidFrames2_1,4); % cam2
numf3 = size(vidFrames3_1,4); % cam3



%% get the video data

for k = 1:numf1 % get data for cam1
    mov1(k).cdata = vidFrames1_1(:,:,:,k); mov1(k).colormap = [];
end
for k = 1:numf2 % get data for cam2
    mov2(k).cdata = vidFrames2_1(:,:,:,k); mov2(k).colormap = [];
end
for k = 1:numf3 % get data for cam3
     mov3(k).cdata = vidFrames3_1(:,:,:,k); mov3(k).colormap = [];
end  

X1 = []; X2 = []; X3 = []; Y1 = []; Y2 = []; Y3 = [];

% RGB to grayscale
% find points of maximum intensity
for k = 1:numf1
    X=frame2im(mov1(k));
    temp = rgb2gray(X);
    Vid1org(:,:,k) = temp;
    temp(:,1:320) = 0; temp(:,380:end) = 0;
    temp(1:200,:) = 0;
    Vid1(:,:,k) =  temp;
    [MAX idx] = max(temp(:));
    [x1 y1] = ind2sub(size(temp), idx);
    X1 = [X1 x1]; Y1 = [Y1 y1];
end

for k = 1:numf2
    X=frame2im(mov2(k));
    temp = rgb2gray(X);
    Vid2org(:,:,k) = temp;
    temp(:,1:260) = 0; temp(:,330:end) = 0;
    Vid2(:,:,k) =  temp;
    [MAX idx] = max(temp(:));
    [x2 y2] = ind2sub(size(temp), idx);
    X2 = [X2 x2]; Y2 = [Y2 y2];
end 

for k = 1:numf3
    X=frame2im(mov3(k));
    temp = rgb2gray(X);
    Vid3org(:,:,k) = temp;
    temp(1:250,:) = 0; temp(310:end,:) = 0;
    temp(:, 1:260) = 0;
    Vid3(:,:,k) =  temp; 
    [MAX idx] = max(temp(:));
    [x3 y3] = ind2sub(size(temp), idx);
    X3 = [X3 x3]; Y3 = [Y3 y3];
end


%% plot X and Y position data
figure (1)
subplot(3,2,1), plot(Y1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]); ylabel('Position');
title('Camera 1 - X Position');
subplot(3,2,2), plot(X1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); 
title('Camera 1 - Y Position');
subplot(3,2,3), plot(X2, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); ylabel('Position');
title('Camera 2 - X Position');
subplot(3,2,4), plot(Y2, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]); 
title('Camera 2 - Y Position');
subplot(3,2,5), plot(X3, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); ylabel('Position'); xlabel('Frame');
title('Camera 3 - X Position');
subplot(3,2,6), plot(Y3, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]);  xlabel('Frame');
title('Camera 3 - Y Position');
hold off


%% put all data in matrix to do SVD
A = [X1(1:200); Y1(1:200); X2(1:200); 
     Y2(1:200); X3(1:200); Y3(1:200)];

[u s v] = svd(A);

Asvd = u(:,1:2)*s(1:2,1:2)*v(:,1:2)';

Xsvd1 = Asvd(1,:); Ysvd1 = Asvd(2,:);
Xsvd2 = Asvd(3,:); Ysvd2 = Asvd(4,:);
Xsvd3 = Asvd(5,:); Ysvd3 = Asvd(6,:);


%% plot two-mode recreation
figure (2)
subplot(3,2,1), plot(Ysvd1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]); ylabel('Position');
title('Camera 1 - X Position');
subplot(3,2,2), plot(Xsvd1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); 
title('Camera 1 - Y Position');
subplot(3,2,3), plot(Xs2, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); ylabel('Position');
title('Camera 2 - X Position');
subplot(3,2,4), plot(Ysvd2, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]); 
title('Camera 2 - Y Position');
subplot(3,2,5), plot(Xsvd3, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); ylabel('Position'); xlabel('Frame');
title('Camera 3 - X Position');
subplot(3,2,6), plot(Ysvd3, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]);  xlabel('Frame');
title('Camera 3 - Y Position');
hold off
%%
figure(3)
plot(diag(s)*100/sum(diag(s)), 'ko--', 'MarkerSize', 10'); hold on;
axis([1 6 0 100]); 
xlabel('Singular Values', 'FontSize',14); 
ylabel('% of Energy', 'FontSize',14);
title('Singular Values (case 1)','FontSize', 18);
