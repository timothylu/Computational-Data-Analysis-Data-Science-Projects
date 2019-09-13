clear all; close all; clc

cam_num = 4; % for convenience of making plots

cropsize11 = [300 160 110 240];
cropsize12 = [230 70 150 310];
cropsize13 = [120 280 130 220];
cropsize21 = [280 220 150 200];
cropsize22 = [180 60 220 340];
cropsize23 = [150 280 150 220];
cropsize31 = [260 230 140 170];
cropsize32 = [220 150 280 300];
cropsize33 = [150 280 170 220];
cropsize41 = [300 250 90 200];
cropsize42 = [200 100 270 260];
cropsize43 = [180 320 170 200];

% load data
if cam_num == 1
    load('cam1_1.mat');
    load('cam2_1.mat');
    load('cam3_1.mat');
    v1 = vidFrames1_1; v2 = vidFrames2_1; v3 = vidFrames3_1;
    clear vidFrames1_1 vidFrames2_1 vidFrames3_1;
    cropsize1 = cropsize11; cropsize2 = cropsize12; cropsize3 = cropsize13;
elseif cam_num == 2
    load('cam1_2.mat');
    load('cam2_2.mat');
    load('cam3_2.mat');
    v1 = vidFrames1_2; v2 = vidFrames2_2; v3 = vidFrames3_2;
    clear vidFrames1_2 vidFrames2_2 vidFrames3_2;
    cropsize1 = cropsize21; cropsize2 = cropsize22; cropsize3 = cropsize23;
elseif cam_num == 3
    load('cam1_3.mat');
    load('cam2_3.mat');
    load('cam3_3.mat');
    v1 = vidFrames1_3; v2 = vidFrames2_3; v3 = vidFrames3_3;
    clear vidFrames1_3 vidFrames2_3 vidFrames3_3;
    cropsize1 = cropsize31; cropsize2 = cropsize32; cropsize3 = cropsize33;
else
    load('cam1_4.mat');
    load('cam2_4.mat');
    load('cam3_4.mat');
    v1 = vidFrames1_4; v2 = vidFrames2_4; v3 = vidFrames3_4;
    clear vidFrames1_4 vidFrames2_4 vidFrames3_4;
    cropsize1 = cropsize41; cropsize2 = cropsize42; cropsize3 = cropsize43;
end



%% Analysis on frame data

% number of frames
nf1 = size(v1,4);

% crop video and convert to RGB
for j=nf1:-1:1
    cropped(:,:,:,j) = imcrop(v1(:,:,:,j), cropsize1);
    grayscale(:,:,:,j) = rgb2gray(cropped(:,:,:,j));
end

Rave1 = zeros(nf1);
Cave1 = zeros(nf1);

% track the can
for j = 1:nf1
    j_frame = double(grayscale(:,:,j));
    [row_1,col_1] = find(j_frame >= 245);
    Rave1(j)=mean(row_1);
    Cave1(j)=mean(col_1);
end

clear cropped grayscale;

% number of frames
nf2 = size(v2,4);
for j=nf2:-1:1
    cropped(:,:,:,j) = imcrop(v2(:,:,:,j), cropsize2);
    grayscale(:,:,:,j) = rgb2gray(cropped(:,:,:,j));
end

Rave2 = zeros(nf2);
Cave2 = zeros(nf2);

% track the can
for j = 1:nf2
    j_frame = double(grayscale(:,:,j));
    [row_2,col_2] = find(j_frame >= 245);
    Rave2(j)=mean(row_2);
    Cave2(j)=mean(col_2);
end

clear cropped grayscale;

% number of frames
nf3 = size(v3,4);
for j=nf3:-1:1
    rot_vid(:,:,:,j) = imrotate(v3(:,:,:,j), -90);
    cropped(:,:,:,j) = imcrop(rot_vid(:,:,:,j), cropsize3);
    grayscale(:,:,:,j) = rgb2gray(cropped(:,:,:,j));
end


%%
Rave3 = zeros(nf3);
Cave3 = zeros(nf3);

% track the can
for j = 1:nf3
    j_frame = double(grayscale(:,:,j));
    [row_3,col_3] = find(j_frame >= 245);
    Rave3(j)=mean(row_3);
    Cave3(j)=mean(col_3);
end

clear cropped grayscale;

%% plot positions

axi = [0 200 0 300];

figure(1);
subplot(3,1,1)
min_frames1=min(size(Rave1,2));
t1=1:min_frames1;
plot(t1, Rave1(1:min_frames1), 'k', t1 ,Cave1(1:min_frames1), 'k--')
axis(axi); 
title('Position vs. Time Case 4', 'FontSize', 18)
legend('Y','X');

subplot(3,1,2)
min_frames2=min(size(Rave2,2));
t2=1:min_frames2;
plot(t2,Rave2(1:min_frames2),'k',t2,Cave2(1:min_frames2),'k--')
axis(axi); 
ylabel('Position', 'FontSize', 14)

subplot(3,1,3)
min_frames3=min(size(Rave3,2));
t3=1:min_frames3;
plot(t3,Rave3(1:min_frames3),'k',t3,Cave3(1:min_frames3),'k--')
axis(axi); 
xlabel('Frame', 'FontSize', 14)

%% SVD
data=[Rave1(1:200); Cave1(1:200); Rave2(1:200);
Cave2(1:200); Rave3(1:200); Cave3(1:200)];


% handle NaN values
data(isnan(data))=0;
[u, s, v] = svd(data, 'econ');
sigma = diag(s);
energy = sigma/sum(sigma);


%% plot singular values
figure(2)
subplot(2,1,1)
plot(energy, 'ko','LineWidth', [1.4])
title('Singular Values Case 4', 'FontSize', 18)
ylabel('Variance (%)', 'FontSize', 14)
subplot(2,1,2)

semilogy(energy, 'ko','LineWidth', [1.4])
ylabel('Variance (log)', 'FontSize', 14)
xlabel('Singular Values', 'FontSize', 14)

