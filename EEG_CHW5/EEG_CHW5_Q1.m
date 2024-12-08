close all
clear
% Part 1
load('ElecPosXYZ') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ; % Radius of diffetent layers
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3]; 
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

% Calculating Gain Matrix
Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
R = 9.2; % Radius of head

% Part 2
% Plotting dipoles and electrodes
figure
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), [], '.')
hold on
for i=1:length(ElecPos)
    electrode = ElecPos{1, i};
    scatter3(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),'*', 'g')
    text(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),electrode.Name)
end
hold off
title("Dipoles and Electrodes Positions.png")
saveas(gcf, 'dipoles and electrodes.png')

% Part 3

N = size(LocMat, 2);
forward_inverse(1200, 'Surface', LocMat, GainMat, ElecPos);

figure
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), [], '.')
hold on
for i=1:length(ElecPos)
    electrode = ElecPos{1, i};
    scatter3(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),'*', 'g')
    text(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),electrode.Name)
end
hold off
title("Dipoles and Electrodes Positions.png")
saveas(gcf, 'dipoles and electrodes.png')


forward_inverse(500, 'Deep', LocMat, GainMat, ElecPos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch Dipoles
%%
figure
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), [], '.')
hold on
for i=1:length(ElecPos)
    electrode = ElecPos{1, i};
    scatter3(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),'*', 'g')
    text(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),electrode.Name)
end
hold off
name = "Patch";
patch_indexes = find_nearest(LocMat(:, 1100), LocMat, 2);
n_dipoles = length(patch_indexes);

X = LocMat(:, patch_indexes);
e = X/norm(X, 1);
%e = reshape(e, [], 1);
hold on
scatter3(X(1, :), X(2, :), X(3, :), 'red', 'filled', 'o')
for i=1:size(X,2)
    quiver3(X(1, i), X(2, i), X(3, i), X(1,i)+e(1,i)/10, X(2,i)+e(2,i)/10, X(3,i)+e(3,i)/10)
end
hold off
title("Patch"+" Dipole.png")
saveas(gcf, ["Patch" + 'dipole.png'])

% Part 4
% Loading spikey signal
Interictal = double(load('Interictal.mat').Interictal);
q =  Interictal(1, :);
%Q = e * Interictal(1, :);
% Calculating potentials
M = 0;
for i=1:length(patch_indexes)
    dp_index = patch_indexes(i);
    M = M + GainMat(:, 3*dp_index-2:3*dp_index) * e(:, i) * Interictal(1, :);
end


fs = 256;
%Plotting potentials
%disp_eeg(M, [], fs)
%title("Electrode Potentials " + "Patch")
%saveas(gcf, "Patch"+" potentials.png")

% Part 5
T = [1722, 1785, 5385, 5451, 6321, 8328, 9365, 9439]; % seizure moments
Spike_avg = zeros([size(M, 1), 7]);
for t=T
    Spike_avg = Spike_avg + M(:, t-3:t+3);
end
Spike_avg = sum(Spike_avg, 2)/length(T)/7;

R = 9.2;
figure
hold on
color = (Spike_avg - min(Spike_avg(:)));
color = color/max(color);
for i=1:length(ElecPos)
    electrode = ElecPos{1, i};
    scatter3(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),60,color(i), "filled")
    text(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),electrode.Name)
    colorbar
    %colormap jet
end
view(3)
hold off
title(name + " - Electrodes Potentials")
saveas(gcf, name + " - Electrodes Potentials.png")

figure
Display_Potential_3D(R, Spike_avg)
title("3D Potentials "+name)
saveas(gcf, name+" 3D Potentials.png")

% Part 6
% Solving inverse problem by MNE
alpha = 1;
Q_mne = GainMat'* inv(GainMat*GainMat'+ alpha * eye(size(GainMat, 1))) * M;

% WMNE
W = W_MNE(GainMat, LocMat);
Q_wmne = inv(W'*W)* GainMat' * inv(GainMat*inv(W'*W)*GainMat' + alpha * eye(size(GainMat, 1)))*M;

% LORETA
W = W_LORETA(GainMat, LocMat);
Q_loreta = inv(W'*W)* GainMat' * inv(GainMat*inv(W'*W)*GainMat' + alpha * eye(size(GainMat, 1)))*M;
     

% Estimating dipole position
Q_mne = reshape(mean(Q_mne.^2, 2), 3, []);
Q_wmne = reshape(mean(Q_wmne.^2, 2), 3, []);
Q_loreta = reshape(mean(Q_loreta.^2, 2), 3, []);

[B, dipoles_estimated ] = maxk(sum(Q_mne), n_dipoles);
[B, dipoles_estimated ] = maxk(sum(Q_wmne), n_dipoles);
[B, dipoles_estimated ] = maxk(sum(Q_mne), n_dipoles);

hold on
scatter3(X(1, :), X(2, :), X(3, :), 'red', 'filled', 'o')
scatter3(LocMat(1, dipoles_estimated), LocMat(2, dipoles_estimated), LocMat(3, dipoles_estimated), ...
    'green', 'filled', 'o')

hold off
%scatter3(LocMat(1, i_max)+e_max(1), LocMat(2, i_max)+e_max(2), LocMat(3, i_max)+e_max(3), 'green', 'filled', 'o')
title("Estimated Dipole "+name)
saveas(gcf, name+" Estimated Dipole.png")

ROC("MNE", Q_mne, patch_indexes);
ROC("WMNE", Q_wmne, patch_indexes);
ROC("LORETA", Q_loreta, patch_indexes);
    


function [error_position, error_direction] = forward_inverse(dipole_number, name, LocMat, GainMat, ElecPos)
    dp_index = dipole_number;
    X0 = LocMat(:, dp_index);
    e0 = X0/norm(X0);
    hold on
    scatter3(X0(1), X0(2), X0(3), 'red', 'filled', 'o')
    quiver3(X0(1), X0(2), X0(3),X0(1)+e0(1)/10, X0(2)+e0(2)/10, X0(3)+e0(3)/10)
    hold off
    title(name+" Dipole.png")
    saveas(gcf, [name ' dipole.png'])
    
    % Part 4
    % Loading spikey signal
    Interictal = double(load('Interictal.mat').Interictal);
    Q = e0 * Interictal(1, :);
    % Calculating potentials
    M = GainMat(:, 3*dp_index-2:3*dp_index) * Q;
    fs = 256;
    %Plotting potentials
    %disp_eeg(M, [], fs)
    %title("Electrode Potentials " + name)
    %saveas(gcf, name+" potentials.png")
    
    % Part 5
    T = [1722, 1785, 5385, 5451, 6321, 8328, 9365, 9439]; % seizure moments
    Spike_avg = zeros([size(M, 1), 7]);
    for t=T
        Spike_avg = Spike_avg + M(:, t-3:t+3);
    end
    Spike_avg = sum(Spike_avg, 2)/length(T)/7;
    
    R = 9.2;
    figure
    hold on
    color = (Spike_avg - min(Spike_avg(:)));
    color = color/max(color);
    for i=1:length(ElecPos)
        electrode = ElecPos{1, i};
        scatter3(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),60,color(i), "filled")
        text(R*electrode.XYZ(1), R*electrode.XYZ(2), R*electrode.XYZ(3),electrode.Name)
        colorbar
        %colormap jet
    end
    view(3)
    hold off
    title(name + " - Electrodes Potentials")
    saveas(gcf, name + " - Electrodes Potentials.png")
    
    figure

    Display_Potential_3D(R, Spike_avg)
    hold on
    scatter3(X0(1), X0(2), X0(3), 'red', 'filled', 'o')
    %quiver3(X(1), X(2), X(3),X(1)+e(1), X(2)+e(2), X(3)+e(3))
    title("3D Potentials "+name)
    hold off
    saveas(gcf, name+" 3D Potentials.png")
    % Part 6
    % Solving inverse problem by MNE
    alpha = 1;
    Q_mne = GainMat'* inv(GainMat*GainMat'+ alpha * eye(size(GainMat, 1))) * M;

    % WMNE
    W = W_MNE(GainMat, LocMat);
    Q_wmne = inv(W'*W)* GainMat' * inv(GainMat*inv(W'*W)*GainMat' + alpha * eye(size(GainMat, 1)))*M;

    % LORETA
    W = W_LORETA(GainMat, LocMat);
    Q_loreta = inv(W'*W)* GainMat' * inv(GainMat*inv(W'*W)*GainMat' + alpha * eye(size(GainMat, 1)))*M;
     
    % Part 7
    
    % Estimating dipole position
    i_mne = find_largest_dipole(Q_mne);
    i_wmne = find_largest_dipole(Q_wmne);
    i_loreta = find_largest_dipole(Q_loreta);
    
    % Part 8
    error_position = @(x, x_pred) norm(x-x_pred);
    error_direction = @(e, e_pred) acos(e'*e_pred)*180/pi;
    
    X = LocMat(:, i_mne);
    e = X/norm(X);
    fprintf("Position Error mne = %d cm - %s\n", error_position(X0, X), name)
    fprintf("Direction Error mne = %d degrees %s\n", error_direction(e0, e), name)
    X = LocMat(:, i_wmne);
    e = X/norm(X);
    fprintf("Position Error wmne = %d cm - %s\n", error_position(X0, X), name)
    fprintf("Direction Error wmne = %d degrees %s\n", error_direction(e0, e), name)
    X = LocMat(:, i_loreta);
    e = X/norm(X);
    fprintf("Position Error loreta = %d cm - %s\n", error_position(X0, X), name)
    fprintf("Direction Error loreta = %d degrees %s\n", error_direction(e0, e), name)
    
    % Part 9 
    
    % Define the objective function
    dir = @(phi, theta) [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]';

    m = mean(abs(M), 2);

    fun = @(q) -norm( q(4)*GainMat(:, 3*find_index(q(1:3), LocMat)-2: 3*find_index(q(1:3), LocMat)) ...
        *dir(q(5), q(6)) - m );

    % Chromosome: [x y z amplitude phi theta]
    
    % Define the number of variables
    nvars = 6;
    
    % Define the lower and upper bounds for the variables
    lb = [min(LocMat(1,:)), min(LocMat(2,:)), min(LocMat(3,:)), 0, 0, 0];
    ub = [max(LocMat(1,:)), max(LocMat(2,:)), max(LocMat(3,:)), norm(m)*30, 2*pi, pi];
    
    % Set the options for the genetic algorithm
    options = optimoptions("ga","MaxGenerations", 200, "PopulationSize", 100 , ...
        "PlotFcn",["gaplotbestf"]);
    
    % Run the genetic algorithm
    [x,fval] = ga(fun,nvars,[],[],[],[],lb,ub,[],[],options);


    % Display the results
    fprintf("GA Position Error = %d cm - %s\n", error_position(X0, x(1:3)), name)
    fprintf("GA Direction Error = %d degrees %s\n", error_direction(dir(x(5), x(6)), e0), name)
    hold off
end


function index = find_largest_dipole(Q)
    Q = mean(Q.^2, 2);
    Q = reshape(Q, 3, []);
    [B, index] = max(sum(Q));
end


function index = find_index(x, LocMat)
    distance = sqrt((LocMat(1, :)-x(1)).^2 + (LocMat(2, :)-x(2)).^2 + (LocMat(3, :)-x(3)).^2);
    index = find(min(distance)==distance);
end

function indexes = find_nearest(x, LocMat, radius)
    distance = sqrt((LocMat(1, :)-x(1)).^2 + (LocMat(2, :)-x(2)).^2 + (LocMat(3, :)-x(3)).^2);
    indexes = find(distance<radius);
end

function W = W_MNE(GainMat, LocMat)
    omega = zeros(size(LocMat, 2));
    for b=1:size(LocMat, 2)
        omega(b, b) = sqrt(sum(GainMat(:, 3*b-2:3*b).*GainMat(:, 3*b-2:3*b), 'all'));
    end
    W = kron(omega, eye(3));
end

function W = W_LORETA(GainMat, LocMat)
    P = size(LocMat, 2);
    omega = zeros(P);
    for b=1:P
        omega(b, b) = sqrt(sum(GainMat(:, 3*b-2:3*b).*GainMat(:, 3*b-2:3*b), 'all'));
    end
    
    A1 = zeros(P);
    
    d = 1;
    for i=1:P-1
        for j=i+1:P
            dist = norm(LocMat(:, i)-LocMat(:, j));
            if (d-0.1<dist) && (d+0.1>dist)
                A1(i, j) = 1/6;
                A1(j, i) = 1/6;
            end
          
        end
    end
    A0 = inv(diag(A1*ones(P,1)))*A1;
    B = 6/d^2 * (kron(A0, eye(3)) - eye(3*P));
    W = kron(omega, eye(3)) *(B'*B) *kron(omega, eye(3));
end

function AUC = ROC(name, Q, patch_indexes)
    Q = sum(Q);
    TH = (0:0.05:1)*max(Q, [], "all");
    TPR = [];
    FPR = [];
    for th=TH
        index_estimate = find(Q>=th);
        TP = sum(ismember(index_estimate, patch_indexes));
        FP = length(index_estimate) - TP;
    
        TPR(end+1) = TP/length(patch_indexes);
        FPR(end+1) = FP/(size(Q,2)-length(patch_indexes));
    end
    figure 
    plot(FPR, TPR)
    xlabel("FPR")
    ylabel("TPR")
    title(name + " ROC")
    saveas(gcf, name + " ROC.png")
end
