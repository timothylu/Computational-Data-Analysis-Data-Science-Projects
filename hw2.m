%% setup

clear all; close all; clc;

load handel
v = y'/2;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

v = v(1:length(v) - 1);
L = 9;
n = length(v);
t = (1:length(v))/Fs;


k=(2*pi/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
%p8 = audioplayer(v,Fs);
%playblocking(p8);  



%% gabor filtering on Handel
Vgt_spec=[];

step = L/90;
tslide=0:.01:L;

a = 500;

for j=1:length(tslide)
    tf = t-tslide(j);

    g=exp(-a*tf.^2); % Gaussian
    
    Vg=g.*v; Vgt=fft(Vg);
        
    Vgt_spec=[Vgt_spec; abs(fftshift(Vgt))];
    
    subplot(3,1,1), plot(t,v,'k',t,g,'r')
    axis([0 L -0.5 1]);
    subplot(3,1,2), plot(t,Vg,'k')
    axis([0 L -0.5 1]);
    subplot(3,1,3), plot(t,abs(fftshift(Vgt))/max(abs(Vgt)),'k')
    drawnow
    pause(0.01)
end


subplot(1,1,1), pcolor(tslide,ks,Vgt_spec.'), shading interp
colormap(hot)
title("mexican hat filter")
xlabel("time (s)");
ylabel("frequency (Hz)");

%% Mexican Hat

Vgt_spec=[];

step = L/90;
tslide=0:.1:L;

a = 500;
b = 1000;

for j=1:length(tslide)
    tf = t-tslide(j);

    g=(1-b*tf.^2).*exp(-a*tf.^2); % Mexican hat
    
    Vg=g.*v; Vgt=fft(Vg);
        
    Vgt_spec=[Vgt_spec; abs(fftshift(Vgt))];
    
    subplot(3,1,1), plot(t,v,'k',t,g,'r')
    axis([0 L -0.5 1]);
    title('Unfiltered Data')
    xlabel('Amplitude')
    ylabel('Time(s)')
    subplot(3,1,2), plot(t,Vg,'k')
    axis([0 L -0.5 1]);
    title('Filtered Data')
    xlabel('Amplitude')
    ylabel('time')
    subplot(3,1,3), plot(t,abs(fftshift(Vgt))/max(abs(Vgt)),'k')
    title('Normalized FFT')
    xlabel('Frequency')
    ylabel('Amplitude')
    drawnow
    pause(0.01)
end



subplot(1,1,1),pcolor(tslide,ks,Vgt_spec.'), shading interp
colormap(hot)
title("mexican hat filter");

xlabel("time (s)");
ylabel("frequency (Hz)");

%% Shannon

Vgt_spec=[];

step = L/90;
tslide=0:.1:L;

width=0.25;
idxWidth = round(width/2/(L/n));

for j=1:length(tslide)
    tf = t-tslide(j);

    g=abs(tf) <= width/2;
    
    Vg=g.*v; Vgt=fft(Vg);
        
    Vgt_spec=[Vgt_spec; abs(fftshift(Vgt))];
    
    subplot(3,1,1), plot(t,v,'k',t,g,'r')
    axis([0 L -0.5 1]);
    subplot(3,1,2), plot(t,Vg,'k')
    axis([0 L -0.5 1]);
    subplot(3,1,3), plot(t,abs(fftshift(Vgt))/max(abs(Vgt)),'k')
    drawnow
    pause(0.01)
end



subplot(1,1,1),pcolor(tslide,ks,Vgt_spec.'), shading interp
colormap(hot)
title("Shannon step filter");

xlabel("time (s)");
ylabel("frequency (Hz)");


%% part 2

clear all; close all; clc;

tr_piano=16; % record time in seconds
y=audioread('music1.wav'); Fs=length(y)/tr_piano;
plot((1:length(y))/Fs,y)
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow
%p8 = audioplayer(y,Fs); playblocking(p8);
%figure(2)

v = y'/2;
v = v(1:length(v));
L = length(v)/Fs;
n = length(v);
t = (1:length(v))/Fs;

k=(2*pi/(L))*[0:(n/2-1) (-n/2:-1)]; ks=fftshift(k);

%% gabor filtering 
Vgt_spec=[];

step = .1;
tslide=0:step:L;

tau = 500;


for j=1:length(tslide)
    g=exp(-tau*(t-tslide(j)).^2); % Gaussian
    
    Vg=g.*v; Vgt=fft(Vg);
        
    Vgt_spec=[Vgt_spec; abs(fftshift(Vgt))];
    
    subplot(3,1,1), plot(t,v,'k',t,g,'r')
    axis([0 L -0.5 1]);
    subplot(3,1,2), plot(t,Vg,'k')
    axis([0 L -0.5 1]);
    subplot(3,1,3), plot(ks,abs(fftshift(Vgt)),'k')
    drawnow
    pause(0.01)
end
%% get frequencies
freq = [];

ksHz = ks/(2*pi);

for ii=1:length(Vgt_spec)
    [maximum, index] = max(Vgt_spec(ii));
    if maximum > 0
        freq = [freq; index];
    end
end
freq

%% Spectogram - Gaussian


%Vgt_spec = Vgt_spec./max(abs(Vgt_spec));

subplot(1,1,1), pcolor(tslide,ks/(2*pi),Vgt_spec.'), shading interp
%color(hot)

set(gca,{'Ylim', 'Fontsize'}, {[100 1000],  [14]}) 

title("Mary Had A Little Lamb (Piano) With Overtones")
xlabel("time (s)");
ylabel("frequency (Hz)");
%% Score the music

pcolor(tslide,ks/(2*pi),Vgt_spec.'), shading interp
colormap(hot);
hold on
pb= plot([0 16],[246.94	246.94],'c') %B_3 freq
hold on;
pc = plot([0 16],[261.63 261.63],'b'); %C_4 freq
hold on;
pcs=plot([0 16],[277.18 277.18],'g') %C#_4 freq
hold on
pd=plot([0 16],[293.66	 293.66],'y'); %D_4 freq
hold on;
pds=plot([0 16],[311.13	 311.13	],'r') %D#4 freq
hold on;
pe = plot([0 16],[329.63 329.63],'w') %E_4 freq

set(gca,'Ylim', [200 350],'Fontsize',[14]) 
title("Score of Mary Had A Little Lamb (Piano)")
xlabel("time (s)");
ylabel("frequency (Hz)");

lgd = legend([pe, pds, pd, pcs, pc, pb],'E4','D#4','D4','C#4','C4','B3');
lgd.FontSize = 10;
lgd.Title.String='Note Frequencies';

%% Recorder
clear all; close all; clc

tr_rec=14; % record time in seconds
y=audioread('music2.wav'); Fs=length(y)/tr_rec;
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
p8 = audioplayer(y,Fs); playblocking(p8);

v = y'/2;
v = v(1:length(v));
L = length(v)/Fs;
n = length(v);
t = (1:length(v))/Fs;

k=(2*pi/(L))*[0:(n/2-1) (-n/2:-1)]; ks=fftshift(k);

%% gabor filtering 
Vgt_spec=[];

step = .2;
tslide=0:step:L;

tau = 500;


for j=1:length(tslide)
    g=exp(-tau*(t-tslide(j)).^2); % Gaussian
    
    Vg=g.*v; Vgt=fft(Vg);
        
    Vgt_spec=[Vgt_spec; abs(fftshift(Vgt))];
    
    subplot(3,1,1), plot(t,v,'k',t,g,'r')
    axis([0 L -0.5 1]);
    subplot(3,1,2), plot(t,Vg,'k')
    axis([0 L -0.5 1]);
    subplot(3,1,3), plot(ks,abs(fftshift(Vgt)),'k')
    drawnow
    pause(0.01)
end


%% Spectogram - Gaussian


subplot(1,1,1), pcolor(tslide,ks/(2*pi),Vgt_spec.'), shading interp
set(gca,{'Ylim', 'Fontsize'}, {[500 3500],  [14]}) 

title("Mary Had A Little Lamb (Recorder) With Overtones")
xlabel("time (s)");
ylabel("frequency (Hz)");

%% Score the music
pcolor(tslide,ks/(2*pi),Vgt_spec.'), shading interp
colormap(hot);
hold on;
pb = plot([0 16],[783.99 783.99],'c') %G4 freq
hold on;
pc= plot([0 16],[830.61	830.61],'b') %G#4 freq
hold on;
pcs = plot([0 16],[880.00 880.00],'g'); %A4 freq
hold on;
pd=plot([0 16],[932.33 932.33],'y') %A#4 freq
hold on
pds=plot([0 16],[987.77 987.77],'r'); %B4 freq
hold on;
pe=plot([0 16],[1046.50	1046.50],'w') %C5 freq

set(gca,'Ylim', [700 1100],'Fontsize',[14]) 
title("Score of Mary Had A Little Lamb (Recorder)")
xlabel("time (s)");
ylabel("frequency (Hz)");

lgd = legend([pe,pds,pd,pcs,pc,pb],'C6','B5','A#5','A5','G#5','G5');
lgd.FontSize = 10;
lgd.Title.String='Note Frequencies';
