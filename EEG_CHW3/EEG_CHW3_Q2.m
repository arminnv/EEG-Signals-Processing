close all
% Loading data
data = load('Ex2.mat');
X_org = data.X_org;
fs = 250;
N = size(X_org,2);
P = size(X_org, 1);
t = (0:size(X_org, 2)-1)/fs;

% Displaying signal and noise
disp_eeg(X_org, 4, fs, Electrodes.labels, "Original Signal")
saveas(gcf, "original signal.png")

T_on = zeros(1, N);
threshold = 1.5;
for i=1:N
    if abs(X_org(13, i))>threshold
        T_on(i)=1;
    else
        T_on(i)=0;
    end
end

SNRS = [-10, -20];
noise_4 = data.X_noise_4;
noise_5 = data.X_noise_5;
X_noises={noise_4, noise_5};
ind = 3;
for nx=1:2
    ind = ind+1;
    X_noise = X_noises{nx};
    for k=1:length(SNRS)
        SNR = SNRS(k);
        % Calculating sigma
        sigma_2 = sumsqr(X_org)/sumsqr(X_noise) * 10^(-SNR/10); 
        % Adding noise to signal
        X = X_org + X_noise * sqrt(sigma_2);
    
        % GEVD
        Cs = zeros(P);
        Cx = zeros(P);
    
        for i=1:N
            Cx = Cx + X(:, i)*X(:, i)'/N;
        end
        
        for i=1:N
            if T_on(i)>0
                Cs = Cs + X(:, i)*X(:, i)'/N;
            end
        end
         
        [W, D] = eig(Cs, Cx);
        [D,I] = sort(abs(diag(D)), 'descend');
        W_GEVD = W(:, I);
        
        y_gevd = W_GEVD'*X;
        y_gevd(2:end, :) = 0;
        X_den_gevd = inv(W_GEVD')*y_gevd;
        
        % DSS
        [U, S, V] = svd(X*X');
        Z = S^(-1/2) * U' * X;
    
        w_dss = randn(P,1);
    
        for iteration=1:100
            r = w_dss'*Z;
            r_new = r.*T_on;
            w_new = Z*r_new';
            w_dss = w_new/norm(w_new);
        end
        
        X_den_dss = U * S^(1/2) * w_dss * r;
    
        % Calculating error
        RRMSE = sqrt(sumsqr(X_den_gevd - X_org))/sqrt(sumsqr(X_org));
        fprintf('noise %d ,SNR = %d : GEVD RRMSE = %d \n', ind, SNR, RRMSE)
    
        RRMSE = sqrt(sumsqr(X_den_dss - X_org))/sqrt(sumsqr(X_org));
        fprintf('noise %d ,SNR = %d : DSS RRMSE = %d \n', ind, SNR, RRMSE)
        
        % Plotting channel 13 and 24 from each step
        for channel = [13 24]
            figure;
            subplot(3, 1, 1)
            plot(t, X_den_gevd(channel, :))
            title("Denoised")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            subplot(3, 1, 2)
            plot(t, X_org(channel, :))
            title("Original")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            subplot(3, 1, 3)
            plot(t, X(channel, :))
            title("Noisy")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            sgtitle("GEVD - channel " + num2str(channel)+ ", SNR="+num2str(SNR)+ ", noise "+num2str(ind))
            saveas(gcf, "GEVD - channel" + num2str(channel)+ ", SNR="+num2str(SNR) + ", noise "+num2str(ind)+".png")
        end
    
        for channel = [13 24]
            figure;
            subplot(3, 1, 1)
            plot(t, X_den_dss(channel, :))
            title("Denoised")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            subplot(3, 1, 2)
            plot(t, X_org(channel, :))
            title("Original")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            subplot(3, 1, 3)
            plot(t, X(channel, :))
            title("Noisy")
            xlabel("Time(s)")
            ylabel("Amplitude(uV)")
            sgtitle("DSS - channel " + num2str(channel)+ ", SNR="+num2str(SNR) + ", noise "+num2str(ind))
            saveas(gcf, "DSS - channel" + num2str(channel)+ ", SNR="+num2str(SNR) + ", noise "+num2str(ind)+".png")
        end   
    end
end