close all;
clear;
EEG = load('NewEEGSignal.mat').NewEEGSignal;
N = length(EEG); 
fs = 256; 
t = (0:N-1)/fs; 
f = (0:N-1)*fs/N; 
X = abs(fft(EEG));

figure;
subplot(3,1,1);
plot(t,EEG); 
xlabel('Time (s)'); ylabel('Amplitude');
title('Time Domain Signal');
subplot(3,1,2);
plot(f(1:N/2), X(1:N/2));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal');
subplot(3,1,3);
spectrogram(EEG, 'yaxis', [], [], [], fs);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Short-Time Fourier Transform');
saveas(gcf, 'part1.png')

figure;
plot(f(1:N/2+1),abs(X(1:N/2+1)));
xlim([0 80]);
xlabel('Frequency (Hz)')
title('Frequency Domain Signal - limited');
saveas(gcf, 'fft limited.png')

fc = 64;
EEG_filtered = lowpass(EEG,fc, fs);
M = 2;
EEG_ds = downsample(EEG_filtered, M); 

Nds = length(EEG_ds);  
fsds = fs/M; 
tds = (0:Nds-1)/fsds;
fds = (0:Nds-1)*fsds/Nds; 
Xds = abs(fft(EEG_ds)); 

figure;
subplot(3,1,1);
plot(tds,EEG_ds); 
xlabel('Time (s)'); ylabel('Amplitude');
title('Time Domain Signal with Downsampled EEG');
subplot(3,1,2);
plot(fds(1:Nds/2+1),abs(Xds(1:Nds/2+1))); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(3,1,3);
spectrogram(EEG_ds, 'yaxis', [], [], [], fsds); 
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Short-Time Fourier Transform with Downsampled EEG');
saveas(gcf, 'EEG ds.png')

N = length(EEG_ds);
w8 = EEG_ds(1:N/8);
w4 = EEG_ds(1:N/4);
w2 = EEG_ds(1:N/2);
X8 = fft(w8, N/8);
X4 = fft(w4, N/4);
X2 = fft(w2, N/2);
f8 = (0:N/8-1)*fsds/(N/8);
f4 = (0:N/4-1)*fsds/(N/4);
f2 = (0:N/2-1)*fsds/(N/2);

figure;
subplot(5,1,1);
plot(f(1:length(EEG)/4),X(1:length(EEG)/4)); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,2);
plot(fds(1:Nds/2+1),abs(Xds(1:Nds/2+1))); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,3);
plot(f2(1:Nds/4+1),abs(X2(1:Nds/4+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,4);
plot(f4(1:Nds/8+1),abs(X4(1:Nds/8+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,5);
plot(f8(1:Nds/16+1),abs(X8(1:Nds/16+1))); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
sgtitle('DFT without padding')
saveas(gcf, 'dft.png')

figure;
subplot(5,1,1);
plot(f(1:length(EEG)/4),X(1:length(EEG)/4)); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,2);
plot(fds(1:Nds/2+1),abs(Xds(1:Nds/2+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,3);
X2 = fft(w2, N);
plot(fds(1:Nds/2+1),abs(X2(1:Nds/2+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,4);
X4 = fft(w4, N);
plot(fds(1:Nds/2+1),abs(X4(1:Nds/2+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
subplot(5,1,5);
X8 = fft(w8, N);
plot(fds(1:Nds/2+1),abs(X8(1:Nds/2+1))); 
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Domain Signal with Downsampled EEG');
title('Frequency Domain Signal with Downsampled EEG');
sgtitle('DFT with padding')
saveas(gcf, 'dft padding.png')
