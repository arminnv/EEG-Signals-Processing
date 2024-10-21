lose all;
clear;
% Part 1: noise removal of simulated signals
% Loading data
X_noise = load("X_noise.mat").X_noise;
X_org = load("X_org.mat").X_org;
Electrodes = load("Electrodes.mat").Electrodes;
fs = 250;
t = (0:size(X_noise, 2)-1)/fs;

% Displaying signal and noise
disp_eeg(X_org, 4, fs, Electrodes.labels, "Original Signal")
saveas(gcf, "original signal.png")
disp_eeg(X_noise, 4, fs, Electrodes.labels, "Noise Signal")
saveas(gcf, "noise signal.png")

% Target sources
hold_channels = {[2 5 12 26], [15 18 20 30]};
SNRS = [-5 -15];

for i=1:2
    SNR = SNRS(i);
    % Calculating sigma
    sigma_2 = sumsqr(X_org)/sumsqr(X_noise) * 10^(-SNR/10); 
    % Adding noise to signal
    X = X_org + X_noise * sqrt(sigma_2);
    disp_eeg(X, [], fs, Electrodes.labels, "Noisy Signal - SNR="+num2str(SNR))
    saveas(gcf, "noisy signal" + num2str(SNR) +".png")
    
    % Getting independent sources
    [F, W, K] = COM2R(X, 32);
    Z = W*X;
    % Displaying sources
    disp_eeg(Z, [], fs, [], "Sources - SNR="+num2str(SNR))
    saveas(gcf, "Sources - SNR="+num2str(SNR) +".png")

    % Reconstructing signals
    X_den = F(:, hold_channels{i}) * Z(hold_channels{i}, :);
    disp_eeg(X_den, [], fs, Electrodes.labels, "Denoised Signal - SNR="+num2str(SNR))
    saveas(gcf, "Denoised Signal - SNR="+num2str(SNR) +".png")

    % Calculating error
    RRMSE = sqrt(sumsqr(X_den - X_org))/sqrt(sumsqr(X_org));
    fprintf('SNR = %d : RRMSE = %d \n', SNR, RRMSE)
    
    % Plotting channel 13 and 24 from each step
    for channel = [13 24]
        figure;
        subplot(3, 1, 1)
        plot(t, X_den(channel, :))
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
        sgtitle("channel " + num2str(channel)+ ", SNR="+num2str(SNR))
        saveas(gcf, "channel" + num2str(channel)+ ", SNR="+num2str(SNR) +".png")
    end
end
