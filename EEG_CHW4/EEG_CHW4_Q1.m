% EEG CHW4 Q4
clear
close all

% Part 1
eeg = load("ERP_EEG.mat").ERP_EEG;
fs = 240;
L = size(eeg, 1);
Ns = size(eeg, 2);
t = (0:L-1)/fs;

figure
hold on
labels = [];
for n=[100:100:2500]
    % Calculating average signals
    eeg_sync = mean(eeg(:, 1:n), 2);
    plot(t, eeg_sync)
    labels = [labels; n];
end
legend(num2str(labels), Location="bestoutside")
xlabel('Time(s)')
ylabel('ERP')
hold off
title("Syncronized Signals")
saveas(gcf, 'Syncronized Signals.png')

% Part 2
max_abs = zeros(Ns, 1);

eeg_sync = zeros(Ns, 240);
for n=1:Ns
    % Calculating average signals
    eeg_sync(n, :) = mean(eeg(:, 1:n), 2);
    % Finding maximum absolute value of synchtonized signal
    max_abs(n) = max(abs(eeg_sync(n, :)));
end

% Plotting values per number of samples
figure
plot(max_abs)
xlabel('n')
ylabel('Max absolute value')

title("Max Absolute Value - N")
saveas(gcf, 'Max Absolute Value - N.png')

% Part 3
RMS = zeros(Ns-1);
for n=2:Ns
    % RMS error of n and n-1 signals
    RMS(n) = rms(eeg_sync(n) - eeg_sync(n-1));
end

figure
plot(2:Ns, RMS)
xlabel('n')
ylabel('RMS')
title('RMS error')
saveas(gcf, 'rms_error.png');

% Part 4
% Based on the results from part 1-3 best number of samples is 600
N0 = 600;

% Part 5
figure 
hold on
plot(t, mean(eeg(:, 1:N0), 2))
plot(t, mean(eeg, 2))
plot(t, mean(eeg(:, 1:N0/3), 2))
plot(t, mean(eeg(:, randi(Ns, N0)), 2))
plot(t, mean(eeg(:, randi(Ns, N0/3)), 2))
hold off

legend(["N0", "2550", "N0/3", "N0 random", "N0/3 random"])
xlabel('Time(s)')
ylabel('EEG')
title('Part 5')
saveas(gcf, 'Part 5.png');

