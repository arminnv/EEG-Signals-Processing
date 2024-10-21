% EEG CHW4 Q4
close all
%clear

% Part 1
% Loading data
fs = 256;
load("Ex3.mat");
load("AllElectrodes.mat")

ind0 = find(TrainLabel==0);
ind1 = find(TrainLabel==1);

C0 = zeros(30);
C1 = zeros(30);

% Calculating covaricance matrices
for i=ind0
    C0 = C0 + TrainData(:, :, i)*TrainData(:, :, i)'/length(ind0);
end

for i=ind1
    C1 = C1 + TrainData(:, :, i)*TrainData(:, :, i)'/length(ind1);
end

% Getting weights by GEVD
[W, D] = eig(C0, C1);
[D,I] = sort(diag(D), 'descend');
W = W(:, I);

% CSP filters
Wcsp = W;

% Filtering signals
X_filtered = zeros(165, size(Wcsp, 2), 256);

for i=1:165
    X_filtered(i, :, :) = Wcsp'* TrainData(:, :, i);
end

t = (0:size(X_filtered, 3)-1)/fs;
figure

% 9: class 1, 2: class 0, channel C3
ch = 11;

subplot(2, 2, 1)
hold on
plot(t, squeeze(X_filtered(2, end, :)))
plot(t, squeeze(X_filtered(2, 1, :)))
hold off
title("class 0 filtered")
legend(["class 1", "class 0"])
subplot(2, 2, 3)
plot(t, squeeze(TrainData(ch, :, 2)))
title("class 0")

subplot(2, 2, 2)
hold on
plot(t, squeeze(X_filtered(9, end, :)))
plot(t, squeeze(X_filtered(9, 1, :)))
hold off
title("class 1 filtered")
legend(["class 1", "class 0"])
subplot(2, 2, 4)
plot(t, squeeze(TrainData(ch, :, 9)))
title("class 1")

saveas(gcf, "csp filtered signals.png")

% Part 2?
%%
%*******************
electrode_names = ["AFz", "F7", "F3", "Fz", "F4", "F8", "FC3", "FCz", "FC4", ...
    "T7", "C3", "Cz", "C4", "T8", "CP3", "CPz", "CP4", "P7", "P5", "P3", "P1",...
    "Pz", "P2", "P4", "P6", "P8", "PO3", "PO4", "O1", "O2"];
electrodes = [];
for j=1:length(electrode_names)
    for i=1:length(AllElectrodes)
        ind = cell2mat(strfind(electrode_names, AllElectrodes(i).labels));
        if strcmp(electrode_names(j), AllElectrodes(i).labels)
            electrodes = [electrodes, AllElectrodes(i)];
            break
        end
    end
end
figure("WindowState","fullscreen")
subplot(1, 2, 1)
plottopomap([electrodes.X]',[electrodes.Y]',electrode_names', Wcsp(:, 1)/norm(Wcsp(:, 1)));
title("class 0 - substraction")
subplot(1, 2, 2)
plottopomap([electrodes.X]',[electrodes.Y]',electrode_names', Wcsp(:, end)/norm(Wcsp(:, end)));
title("class 1 - feet")

saveas(gcf, "plottomap.png")
%%

% Part 3

n_folds = 4;
number_of_filters = 2*[1, 2, 3, 4, 5, 6, 7, 8, 9];

best_val_accuracy = 0;
m_filter_best = 0;
for m=number_of_filters
    average_accuracy = 0;
    average_accuracy_train = 0;
    for k=1:4
        C0_kv = zeros(30);
        C1_kv = zeros(30);
        N = size(TrainData, 3);
    
        idx_val = 1+(k-1)*int32(N/n_folds):k*int32(N/n_folds);
        idx_train = 1:N;
        idx_train(idx_val) = [];
    
        X_train = TrainData(:, :, idx_train);
        Y_train = TrainLabel(idx_train);
        X_val = TrainData(:, :, idx_val);
        Y_val = TrainLabel(idx_val);
        
        ind0 = find(Y_train==0);
        ind1 = find(Y_train==1);

        % Calculating covaricance matrices
        for i=ind0
            C0_kv = C0_kv + TrainData(:, :, i)*TrainData(:, :, i)'/length(ind0);
        end
        
        for i=ind1
            C1_kv = C1_kv + TrainData(:, :, i)*TrainData(:, :, i)'/length(ind1);
        end
                
        % Getting weights by GEVD
        [W_kv, D] = eig(C0_kv, C1_kv);
        [D,I] = sort(diag(D), 'descend');
        W_kv = W_kv(:, I);
        
        % CSP filters
        Wcsp = [W_kv(:, 1:m/2), W_kv(:, end-m/2+1:end)];
        
        % Filtering signals
        X_filtered = zeros(165, m, 256);
        
        for i=1:165
            X_filtered(i, :, :) = Wcsp'* TrainData(:, :, i);
        end
        
        % Variance of rows as features
        F = var(X_filtered, 1, 3);
        
        % Decision Tree Classifier
        Mdl = fitctree(F(idx_train, :), Y_train);
    
        Y_pred = predict(Mdl,F(idx_val, :))';
        accuracy = mean(Y_pred==Y_val);
        average_accuracy = average_accuracy + accuracy/n_folds;
        %accuracy_train = length(find(predict(Mdl,X_train) == y_train))/length(y_train);
        %average_accuracy_train = average_accuracy_train + accuracy_train/K_folds;
    
    end

    if average_accuracy > best_val_accuracy 
        best_val_accuracy = average_accuracy;
        m_filter_best = m;
    end
end
fprintf("Highest validation accuracy = %d \n", best_val_accuracy)
fprintf("Best number of filters = %d \n", m_filter_best)

% Part 4


% Predicting test labels
Wcsp = [W(:, 1:m_filter_best/2), W(:, end-m_filter_best/2+1:end)];

% Filtering signals
X_train = zeros(165, m_filter_best, 256);

for i=1:165
    X_train(i, :, :) = Wcsp'* TrainData(:, :, i);
end

X_test = zeros(size(TestData, 3), m_filter_best, 256);
for i=1:size(TestData, 3)
    X_test(i, :, :) = Wcsp'* TestData(:, :, i);
end

% Variance of rows as features
F_train = var(X_train, 1, 3);
F_test = var(X_test, 1, 3);

% Decision Tree Classifier
Mdl = fitctree(F_train, TrainLabel);

TestLabel = predict(Mdl, F_test);

