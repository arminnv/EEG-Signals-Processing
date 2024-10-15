close all
% Loading data
data = load('Q1.mat');
X = data.X_org;
fs = 100;
N = size(X,2);
P = size(X, 1);
[U, S, V] = svd(X*X');
Z = S^(-1/2) * U' * X;
plot(Z(1, :))

%%
% Part 1
w1 = randn(P,1);
sgn = 1;
T = 4*fs;
for iteration=1:100
    r1 = w1'*Z;
    r1_new = (r1+circshift(r1, sgn*T))/2;
    sgn = sgn * -1;
    w1_new = Z*r1_new';
    w1 = w1_new/norm(w1_new);
end

X1_den = U * S^(1/2) * w1 * r1;

figure
subplot(3, 1, 1)
plot(time, X1(4, :))
title('X1 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X1_den(4, :))
title('X denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original - channel 4')
xlabel('Time(s)')
sgtitle('DSS - Part 1 - x1')
saveas(gcf, 'dss part1.png')

RRMSE = sqrt(sumsqr(X1_den - X1))/sqrt(sumsqr(X1));
fprintf('Part 1 - RRMSE = %d \n', RRMSE)
%%
% Part 2

t_best = 0;
cor_max = 0;
r1_best = 0;
w1_best = 0;
for t=3:0.1:7
    T = int32(t*fs);
    w1 = randn(P,1);
    sgn = 1;
    for iteration=1:100
        r1 = w1'*Z;
        r1_new = r1+circshift(r1, sgn*T);
        sgn = sgn * -1;
        w1_new = Z*r1_new';
        w1 = w1_new/norm(w1_new);
    end
    cor_coef = corrcoef(r1, circshift(r1, sgn*T));
    if cor_coef(1, 2)>cor_max
        cor_max = cor_coef(1, 2);
        r1_best  = r1;
        t_best = t;
        w1_best = w1;
    end
end

X1_den = U * S^(1/2) * w1_best * r1_best;

figure
subplot(3, 1, 1)
plot(time, X1(4, :))
title('X1 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X1_den(4, :))
title('X denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original - channel 4')
xlabel('Time(s)')
sgtitle('DSS - Part 2 - x1')
saveas(gcf, 'dss part2.png')

RRMSE = sqrt(sumsqr(X1_den - X1))/sqrt(sumsqr(X1));
fprintf('Part 2 - RRMSE = %d \n', RRMSE)
%%
% Part 3

w2 = randn(P,1);

for iteration=1:100
    r2 = w2'*Z;
    r2_new = r2.*T1;
    w2_new = Z*r2_new';
    w2 = w2_new/norm(w2_new);
end

X2_den = U * S^(1/2) * w2 * r2;

figure
subplot(3, 1, 1)
plot(time, X2(4, :))
title('X2 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X2_den(4, :))
title('X denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original')
xlabel('Time(s)')
sgtitle('DSS - Part 3 - x2')
saveas(gcf, 'dss part3.png')

RRMSE = sqrt(sumsqr(X2_den - X2))/sqrt(sumsqr(X2));
fprintf('Part 3 - RRMSE = %d \n', RRMSE)

%%
% Part 4

w2 = randn(P,1);

for iteration=1:100
    r2 = w2'*Z;
    r2_new = r2.*T2;
    w2_new = Z*r2_new';
    w2 = w2_new/norm(w2_new);
end

X2_den = U * S^(1/2) * w2 * r2;

figure
subplot(3, 1, 1)
plot(time, X2(4, :))
title('X2 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X2_den(4, :))
title('X denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original')
xlabel('Time(s)')
sgtitle('DSS - Part 4 - x2')
saveas(gcf, 'dss part4.png')

RRMSE = sqrt(sumsqr(X2_den - X2))/sqrt(sumsqr(X2));
fprintf('Part 4 - RRMSE = %d \n', RRMSE)
%%
% Part 5


w3 = randn(P,1);
for iteration=1:100
    r3 = w3'*Z;
    r3_new = bandpass(r3,[10 15],fs);
    w3_new = Z*r3_new';
    w3 = w3_new/norm(w3_new);
end

X3_den = U * S^(1/2) * w3 * r3;

f = (0:N-1)/fs;

figure
subplot(3, 1, 1)
plot(time, X3(4, :))
title('X3 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X3_den(4, :))
title('X3 denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original')
xlabel('Time(s)')
sgtitle('DSS - Part 5 - x2')
saveas(gcf, 'dss part5.png')

figure
subplot(3, 1, 1)
plot(f, fft(X3(4, :)))
title('X2 - channel 4 (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 2)
plot(f, fft(r3))
title('source (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 3)
plot(f, fft(X(4, :)))
title('X original (Frequency domain)')
xlabel('Frequency(Hz)')
sgtitle('DSS - Part 5 - x3 (Frequency domain)')
saveas(gcf, 'dss part5f.png')

RRMSE = sqrt(sumsqr(X3_den - X3))/sqrt(sumsqr(X3));
fprintf('Part 5 - RRMSE = %d \n', RRMSE)

%%
% Part 6
w3 = randn(P,1);
for iteration=1:100
    r3 = w3'*Z;
    r3_new = bandpass(r3,[5 25],fs);
    w3_new = Z*r3_new';
    w3 = w3_new/norm(w3_new);
end

X3_den = U * S^(1/2) * w3 * r3;

f = (0:N-1)/fs;

figure
subplot(3, 1, 1)
plot(time, X3(4, :))
title('X3 - channel 4')
xlabel('Time(s)')
subplot(3, 1, 2)
plot(time, X3_den(4, :))
title('X3 denoised - channel 4')
xlabel('Time(s)')
subplot(3, 1, 3)
plot(time, X(4, :))
title('X original')
xlabel('Time(s)')
sgtitle('DSS - Part 6 - x2')
saveas(gcf, 'dss part6.png')

figure
subplot(3, 1, 1)
plot(f, fft(X3(4, :)))
title('X2 - channel 4 (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 2)
plot(f, fft(r3))
title('source (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 3)
plot(f, fft(X(4, :)))
title('X original (Frequency domain)')
xlabel('Frequency(Hz)')
sgtitle('DSS - Part 6 - x3 (Frequency domain)')
saveas(gcf, 'dss part5f.png')

RRMSE = sqrt(sumsqr(X3_den - X3))/sqrt(sumsqr(X3));
fprintf('Part 6 - RRMSE = %d \n', RRMSE)
