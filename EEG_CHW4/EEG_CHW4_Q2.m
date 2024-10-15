% EEG CHW4 Q2
clear
close all

% Part A
% Using Power Spectral Density of Signals

% Loading SSVEP data
data = load("SSVEP_EEG.mat");
signal = data.SSVEP_Signal;
events = data.Events;
events_samples = data.Event_samples;

fs = 250;
channels = ["Pz", "Oz", "P7", "P8", "O2", "O1"];

% Bandpass filtering signals
signal = bandpass(signal', [1, 40], fs, "Steepness", 1)';

% Splitting trials
trials = zeros(length(events_samples), size(signal, 1), 5*fs);
for i=1:length(events_samples)
    start = events_samples(i);
    trials(i, :, :) = signal(:, start: start + 5*fs-1);
end

for n=1:15
    % Calculating power spectral density of signals
    [Pxx, w] = pwelch(squeeze(trials(n, : , :))', [], [], [], fs);

    figure
    hold on
    for ch=1:length(channels)
        plot(w, Pxx(:, ch))
    end
    legend(channels)
    xlabel("Frequency(Hz)")
    ylabel("PSD")
    title("trial " + num2str(n))
    hold off
    saveas(gcf, "trial " + num2str(n)+".png")
end

%%
% Part B
% Using Cannonical Correlation Analysis
clear
close all

% Loading SSVEP data
data = load("SSVEP_EEG.mat");
signal = data.SSVEP_Signal;
events = data.Events;
events_samples = data.Event_samples;

fs = 250;
channels = ["Pz", "Oz", "P7", "P8", "O2", "O1"];
f_st = [6.5, 7.35, 8.3, 9.6, 11.61];

% Splitting trials
trials = zeros(length(events_samples), size(signal, 1), 5*fs);
for i=1:length(events_samples)
    start = events_samples(i);
    trials(i, :, :) = signal(:, start: start + 5*fs-1);
end

all_Yf = {};
T = size(trials, 3);
for f_stimuli=f_st
    all_Yf{end+1} = get_Yf(f_stimuli, T, fs);
end

f_pred = zeros(length(events), 1);
for m=1:15
    X = squeeze(trials(m, :, :));
    ro = zeros(length(f_st), 1);
    for n=1:length(f_st)
        [A, B, R] = canoncorr(X', all_Yf{1, n}');
        ro(n) = R(1);
    end
    f_pred(m) = f_st(find(ro==max(ro)));
end

accuracy = mean(f_pred==events');
fprintf("accuracy of prediction = %d \n", accuracy)

% Part B3
% Exclude channels: Pz, P7, P8 - 1, 3, 4
new_channels = [2, 5, 6];

f_pred = zeros(length(events), 1);
for m=1:15
    X = squeeze(trials(m, new_channels, :));
    ro = zeros(length(f_st), 1);
    for n=1:length(f_st)
        [A, B, R] = canoncorr(X', all_Yf{1, n}');
        ro(n) = R(1);
    end
    f_pred(m) = f_st(find(ro==max(ro)));
end

accuracy = mean(f_pred==events');
fprintf("accuracy after removing P channels = %d \n", accuracy)

% Part B3
% Reducing window duration to 2s
new_channels = [2, 5, 6];

f_pred = zeros(length(events), 1);
for m=1:15
    X = squeeze(trials(m, new_channels, :));
    ro = zeros(length(f_st), 1);
    for n=1:length(f_st)
        Yf = all_Yf{1, n};
        [A, B, R] = canoncorr(X(:, 1:T/5*2)', Yf(:, 1:T/5*2)');
        ro(n) = R(1);
    end
    f_pred(m) = f_st(find(ro==max(ro)));
end

accuracy = mean(f_pred==events');
fprintf("accuracy after reducing window duration to 2s = %d \n", accuracy)

function Yf = get_Yf(f_stimuli, T, fs)
    n = round(40/f_stimuli);
    t = (0:T-1)/fs;
    Yf = zeros(2*n, T);
    for i=1:n
        Yf(2*i-1:2*i, :) = [sin(2*pi*i*f_stimuli*t); cos(2*pi*i*f_stimuli*t)];
    end
end
