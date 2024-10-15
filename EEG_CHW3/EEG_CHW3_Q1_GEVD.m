close all
% Loading data
data = load('Q1.mat');
X = data.X_org;
X1 = data.X1;
X2 = data.X2;
X3 = data.X3;
fs = 100;
P = size(X, 1);
N = size(X,2);
time = (0:N-1)/fs;

%%
% Part 1
% GEVD

Cx = zeros(size(X, 1), size(X, 1));

for i=1:N
    Cx = Cx + X(:, i)*X(:, i)'/N;
end

Cs1 = zeros(size(X, 1), size(X, 1));
T = 4*fs;
for i=1:N-T
    Cs1 = Cs1 + X(:, i)*X(:, i+T)'/N;
end
Cs1 = (Cs1 + Cs1')/2;

[W, D1] = eig(Cs1, Cx);

[D1,I] = sort(abs(diag(D1)), 'descend');
disp(D1)
W1 = W(:, I);
w1 = W(:, 1);

y1 = W1'*X;
y1(2:end, :) = 0;
X1_den = inv(W1')*y1;

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
plot(time, y1(1, :))
plot(time, X(4, :))
title('X original - channel 4')
sgtitle('Part 1 - x1')
saveas(gcf, 'gevd part1.png')

RRMSE = sqrt(sumsqr(X1_den - X1))/sqrt(sumsqr(X1));
fprintf('Part 1 - RRMSE = %d \n', RRMSE)

%%
% Part 2
%GEVD

t_best = 3;
w1 = 0;
W1 = 0;
e_max = 0;

for t=3:0.1:7
    T = t*fs;
    Cs1 = zeros(size(X, 1), size(X, 1));
    for i=1:N-T
        Cs1 = Cs1 + X(:, i)*X(:, int32(i+T))'/N;
    end
    Cs1 = (Cs1 + Cs1')/2;
    
    [W, D1] = eig(Cs1, Cx);
    [D1,I] = sort(abs(diag(D1)), 'descend');
    if D1(1)>e_max
        t_best = t;
        e_max = D1(1);
        W1 = W(:, I);
        w1 = W(:, 1);
    end
end

y1 = W1'*X;
y1(2:end, :) = 0;
X1_den = inv(W1')*y1;

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
sgtitle('Part 2 - x1')
saveas(gcf, 'gevd part2.png')

RRMSE = sqrt(sumsqr(X1_den - X1))/sqrt(sumsqr(X1));
fprintf('Part 2 - RRMSE = %d \n', RRMSE)

%%
% Part 3

Cs2 = zeros(size(X, 1), size(X, 1));
T = 4*fs;
L = length(T1);
for i=1:N
    if T1(i)>0
        Cs2 = Cs2 + X(:, i)*X(:, i)'/L;
    end
end
%Cs = (Cs1 + Cs1')/2;

[W, D2] = eig(Cs2, Cx);
[D2,I] = sort(abs(diag(D2)), 'descend');
W2 = W(:, I);
w2 = W(:, 1);

y2 = W2'*X;
y2(2:end, :) = 0;
X2_den = inv(W2')*y2;

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
sgtitle('Part 3 - x2')
saveas(gcf, 'gevd part3.png')

RRMSE = sqrt(sumsqr(X2_den - X2))/sqrt(sumsqr(X2));
fprintf('Part 3 - RRMSE = %d \n', RRMSE)
%%
% Part 4 

Cs2 = zeros(size(X, 1), size(X, 1));
L = length(T2);
for i=1:N
    if T2(i)>0
        Cs2 = Cs2 + X(:, i)*X(:, i)'/L;
    end
end
%Cs = (Cs1 + Cs1')/2;

[W, D2] = eig(Cs2, Cx);
[D2,I] = sort(abs(diag(D2)), 'descend');
W2 = W(:, I);
w2 = W(:, 1);

y2 = W2'*X;
y2(2:end, :) = 0;
X2_den = inv(W2')*y2;

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
sgtitle('Part 4 - x2')
saveas(gcf, 'gevd part4.png')

RRMSE = sqrt(sumsqr(X2_den - X2))/sqrt(sumsqr(X2));
fprintf('Part 4 - RRMSE = %d \n', RRMSE)

%%
% Part 5
Cs3 = zeros(size(X, 1), size(X, 1));
f = (0:N-1)*fs/N; 
Xf = zeros(P, N);
for i=1:P
    Xf(i, :) = fft(X(i, :));
end

Lf = length(Xf);
Sx3 = zeros(P, P);
Sx =  Xf(:, 1)*conj(Xf(:, 1)')/Lf;
figure
plot(f,Xf(1,:))

for i=2:int32(Lf/2)
        Sx = Sx + abs(Xf(:, i)*conj(Xf(:, i)')/Lf + Xf(:, Lf-i+2)*conj(Xf(:, Lf-i+2))')/Lf;
end

for i=int32(10*Lf/fs):int32(15*Lf/fs)
     Sx3 = Sx3 + abs(Xf(:, i)*conj(Xf(:, i)')/Lf + Xf(:, Lf-i+2)*conj(Xf(:, Lf-i+2))')/Lf;
end

[W, D3] = eig(Sx3, Sx);
[D3,I] = sort(abs(diag(D3)), 'descend');
W3 = W(:, I);
w3 = W(:, 1);

y3f = W3'*Xf;
y3f(2:end, :) = 0;
y3 = W3'*X;
y3(2:end, :) = 0;
X3_den = inv(W3')*y3;

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
sgtitle('Part 5 - x2')
saveas(gcf, 'gevd part5.png')

figure
subplot(3, 1, 1)
plot(f, fft(X3(4, :)))
title('X2 - channel 4 (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 2)
plot(f, y3)
title('source (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 3)
plot(f, fft(X(4, :)))
title('X original (Frequency domain)')
xlabel('Frequency(Hz)')
sgtitle('Part 5 - x3 (Frequency domain)')
saveas(gcf, 'gevd part5f.png')

RRMSE = sqrt(sumsqr(X3_den - X3))/sqrt(sumsqr(X3));
fprintf('Part 5 - RRMSE = %d \n', RRMSE)
%%
% Part 6

e_max = 0;
fl_best = 5;
fu_best = 25;
w3 = 0;
W3 =  0;
for fl=5:1:25
    for fu=fl+1:1:25        
        L3 = 0;
        Sx3 = zeros(P, P);
        for i=int32(fl*N/fs):int32(fu*N/fs)
            Sx3 = Sx3 + abs(Xf(:, i)*conj(Xf(:, i)')/Lf + Xf(:, Lf-i+2)*conj(Xf(:, Lf-i+2))')/Lf;
            L3 = L3 + 2;
        end
        Sx3 = Sx3/L3;
        
        [W, D3] = eig(Sx3, Sx);
        [D3,I] = sort(abs(diag(D3)), 'descend');
        if D3(1)>e_max
            fl_best = fl;
            fu_best = fu;
            e_max = D3(1);
            W3 = W(:, I);
            w3 = W(:, 1);         
        end
    end
end

y3f = W3'*Xf;
y3f(2:end, :) = 0;
y3 = W3'*X;
y3(2:end, :) = 0;
X3_den = inv(W3')*y3;

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
sgtitle('Part 6 - x2')
saveas(gcf, 'gevd part6.png')

figure
subplot(3, 1, 1)
plot(f, fft(X3(4, :)))
title('X2 - channel 4 (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 2)
plot(f, y3)
title('source (Frequency domain)')
xlabel('Frequency(Hz)')
subplot(3, 1, 3)
plot(f, fft(X(4, :)))
title('X original (Frequency domain)')
xlabel('Frequency(Hz)')
sgtitle('Part 6 - x3 (Frequency domain)')
saveas(gcf, 'gevd part6f.png')

RRMSE = sqrt(sumsqr(X3_den - X3))/sqrt(sumsqr(X3));
fprintf('Part 6 - RRMSE = %d \n', RRMSE)
