close all; clear all; clc;

% setup/get data from matrix
load('cam1_4.mat');
load('cam2_4.mat');
load('cam3_4.mat');
vidFrames1_1 = vidFrames1_4;
vidFrames2_1 = vidFrames2_4;
vidFrames3_1 = vidFrames3_4;

[m1,n1] = size(vidFrames1_1(:,:,1,1));
[m2,n2] = size(vidFrames2_1(:,:,1,1));
[m3,n3] = size(vidFrames3_1(:,:,1,1));

numFrames1=size(vidFrames1_1,4);
numFrames2=size(vidFrames2_1,4);
numFrames3=size(vidFrames3_1,4);
maxFrames = max([numFrames1 numFrames2 numFrames3]);

for k = 1:maxFrames
    if k <= numFrames1
        mov1(k).cdata = vidFrames1_1(:,:,:,k);
        mov1(k).colormap = [];
        
    end
    if k <= numFrames2
        mov2(k).cdata = vidFrames2_1(:,:,:,k);
        mov2(k).colormap = [];
    end
    if k <= numFrames3
         mov3(k).cdata = vidFrames3_1(:,:,:,k);
         mov3(k).colormap = [];
    end  
end

X1 = []; X2 = []; X3 = [];
Y1 = []; Y2 = []; Y3 = [];
for jj = 1:maxFrames
    
    if jj <= numFrames1
        X=frame2im(mov1(jj));
        a = rgb2gray(X);
        Vid1org(:,:,jj) = a;
        a(:,1:290) = 0; a(:,350:end) = 0;
        a(1:200,:) = 0;
        Vid1(:,:,jj) =  a;
        [Max Ind] = max(a(:));
        [x1 y1] = ind2sub(size(a), Ind);
        X1 = [X1 x1]; Y1 = [Y1 y1];
    end
    
    if jj <= numFrames2
        X=frame2im(mov2(jj));
        a = rgb2gray(X);
        Vid2org(:,:,jj) = a;
        a(:,1:230) = 0; a(:,350:end) = 0;
        Vid2(:,:,jj) =  a;
        [Max Ind] = max(a(:));
        [x2 y2] = ind2sub(size(a), Ind);
        X2 = [X2 x2]; Y2 = [Y2 y2];
    end 
 
    if jj <= numFrames3
        X=frame2im(mov3(jj));
        a = rgb2gray(X);
        Vid3org(:,:,jj) = a;
        a(1:200,:) = 0; a(300:end,:) = 0;
        a(:, 1:260) = 0;
        Vid3(:,:,jj) =  a; 
        [Max Ind] = max(a(:));
        [x3 y3] = ind2sub(size(a), Ind);
        X3 = [X3 x3]; Y3 = [Y3 y3];
    end

end

%% 
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


A = [X1(1:200); Y1(1:200); X2(1:200); 
     Y2(1:200); X3(1:200); Y3(1:200)];

[u s v] = svd(A);
%%
Asvd = u(:,1:2)*s(1:2,1:2)*v(:,1:2)';

Xsvd1 = Asvd(1,:); Ysvd1 = Asvd(2,:);
Xsvd2 = Asvd(3,:); Ysvd2 = Asvd(4,:);
Xsvd3 = Asvd(5,:); Ysvd3 = Asvd(6,:);

%%
figure (2)
subplot(3,2,1), plot(Ysvd1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 640]); ylabel('Position');
title('Camera 1 - X Position');
subplot(3,2,2), plot(Xsvd1, 'k', 'LineWidth', 1.5); 
axis([0 200 0 480]); 
title('Camera 1 - Y Position');
subplot(3,2,3), plot(Xsvd2, 'k', 'LineWidth', 1.5); 
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
title('Singular Values (case 4)','FontSize', 18);




