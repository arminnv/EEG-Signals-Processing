% EEG CHW2 Q1
% Part 1
close all;
clear;
fs = 1000;
t = 0:1/fs:2;
x = chirp(t, 100, 1, 200, 'quadratic');

% Part 2
L = 128; 
w1 = rectwin(L);
w2 = triang(L);
w3 = gausswin(L,5);
w4 = hamming(L,'periodic'); 
%wvtool(w1,w2,w3,w4);
%saveas(gcf, 'wvtool.png')

L = 128; 
Noverlap = 0;
nfft = L; 
fs = 1000; 

figure;
subplot(2,2,1);
spectrogram(x,'yaxis',w1,Noverlap,nfft,fs);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Magnitude');
title('Spectrogram with Rectangular Window');
subplot(2,2,2);
spectrogram(x,'yaxis',w2,Noverlap,nfft,fs);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Magnitude');
title('Spectrogram with Triangular Window');
subplot(2,2,3);
spectrogram(x,'yaxis',w3,Noverlap,nfft,fs);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Magnitude');
title('Spectrogram with Gaussian Window');
subplot(2,2,4);
spectrogram(x,'yaxis',w4,Noverlap,nfft,fs);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Magnitude');
title('Spectrogram with Hamming Window');
saveas(gcf, 'spec-window.png')

Noverlap = [0 64 127];
figure;
for i = 1:length(Noverlap)
    subplot(length(Noverlap),1,i)
    spectrogram(x,w1,Noverlap(i),nfft,fs,'yaxis')
    title(['Spectrogram of x(t) using rectangular window and Noverlap = ' num2str(Noverlap(i))])
end
saveas(gcf, 'rect-nfft.png')

L = [32 128 512]; 
figure;
for i = 1:length(L)
    w = rectwin(L(i));
    Noverlap = L(i)-1; 
    nfft = L(i); 
    subplot(length(L),1,i)
    spectrogram(x,w,Noverlap,nfft,fs,'yaxis')
    title(['Spectrogram of x(t) using rectangular window and L = ' num2str(L(i))])
end
saveas(gcf, 'length.png')

L = 128;
nfft = [L L*2 L*4]; 
figure;
for i = 1:length(nfft)
    subplot(length(nfft),1,i)
    spectrogram(x,w1,L/2,nfft(i),fs,'yaxis')
    title(['Spectrogram of x(t) using rectangular window and nfft = ' num2str(nfft(i))])
end
saveas(gcf, 'nfft.png')


S = myspectrogram(x, 100, 99, 100);
figure;
subplot(2, 1, 1)
image(abs(S));
title('my spectrogram')
subplot(2, 1, 2)
spectrogram(x,rectwin(100),99,100,fs,'yaxis')
title('matlab')
saveas(gcf, 'my spectrogram.png')


function S = myspectrogram(x,L,Noverlap,nfft)
    N = length(x);
    w = rectwin(L);
    Nstep = L - Noverlap;
    
    K = floor((N-Noverlap)/Nstep);
    S = zeros(nfft/2,K);
    
    for k = 1:K
        start = (k-1)*Nstep + 1;
        finish = start + L - 1;
        xw = x(start:finish)'.*w;
        Xk = abs(fft(xw,nfft));
        
        S(:,k) = Xk(1:nfft/2);
    end
    S = flipud(S);
end




