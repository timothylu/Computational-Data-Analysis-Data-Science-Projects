clear all; close all; clc;
load("C:\Users\timot\Desktop\amath 482\Testdata.mat")

%% Set Up

L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%% Averaging the Spectrum

uave = zeros(n,n,n);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    utn=fftn(Un);
    uave=uave+utn;
end

absUave = fftshift(abs(uave));
[maximum, index] = max(absUave(:));

% Center Frequency
[I,J,K] = ind2sub([64,64,64],index)

cx = Kx(I,J,K)
cy = Ky(I,J,K)
cz = Kz(I,J,K)

close all, isosurface(Kx,Ky,Kz,fftshift(abs(uave))./max(abs(uave(:))),0.75)
axis([cx-1 cx+1 cy-1 cy+1 cz-1 cz+1]), grid on, drawnow
xlabel('Kx (Frequency)');
ylabel('Ky (Frequency)');
zlabel('Kz (Frequency)');
title('Averaging of the Spectrum');
ax = gca;
ax.FontSize = 15;


%% Filter
clc

b = 0.1;

filter = fftshift(exp(-b *((Kx-Kx(I, J, K)).^2 + (Ky-Ky(I, J, K)).^2 + (Kz-Kz(I,J,K)).^2)));

xx = zeros(20,1);
yy = zeros(20,1);
zz = zeros(20,1);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    uftn=fftn(Un).*filter;
    unf = ifftn(uftn);
    [maximum,index] = max(abs(unf(:)));
    [I,J,K] = ind2sub([64,64,64],index);
    xx(j) = X(I,J,K);
    yy(j) = Y(I,J,K);
    zz(j) = Z(I,J,K);
end

plot3(xx,yy,zz)
grid on
xlabel('x (position)')
ylabel('y (position)')
zlabel('z (position)')
time = [1:20]; txt = arrayfun(@num2str,time,'UniformOutput',false);
text(xx,yy,zz,txt);
title('Path of the Marble in Fluffy')
ax = gca;
ax.FontSize = 15;

%% part 3

[xx(20), yy(20), zz(20)]

%%
clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);

end

close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Unfiltered Noisy Data');
xlabel('x (position)');
ylabel('y (position)');
zlabel('z (position)');
ax = gca;
ax.FontSize = 15;
