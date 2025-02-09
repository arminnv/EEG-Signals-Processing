close all;
fs = 256;
x = load('NewEEGSignal.mat').NewEEGSignal;
[r, lags] = xcorr(x);  

figure;
plot(lags, r)
xlabel('lag')
ylabel('acorr')
title('autocorrelation of signal')
saveas(gcf, 'autocorrelation.png')

N = length(x); 
r = r(lags>=0);
r(2:end) = 2*r(2:end);
Pxx = real(fft(r))/(N*pi); % power spectral density
Pxx = Pxx(1:N/2);
f = (1:N/2)/N*fs;

figure;
subplot(3, 1, 1)
plot(f, Pxx)
xlabel('frequency(Hz)')
title('our method')
subplot(3, 1, 2)
Xf = periodogram(x);
plot(f, Xf(1:N/2))
xlabel('frequency(Hz)')
title('periodogram')
subplot(3, 1, 3)
plot(pwelch(x))
title('pwelch')
xlabel('frequency(Hz)')
sgtitle('Power Spectral Density')
saveas(gcf, 'PSD.png')
