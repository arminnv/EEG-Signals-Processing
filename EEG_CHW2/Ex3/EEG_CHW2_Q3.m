%% Q2
clc
clear
close all

%Loadind EEG data
load('Electrodes');
load('NewData1');
D1 = EEG_Sig;
load('NewData2');
D2 = EEG_Sig;
load('NewData3');
D3 = EEG_Sig;

fs = 250; %Sampling Freq
ElecName = Electrodes.labels; %Channel labels


%Plotting the Original Data
offset = max(max(abs(D1)))/1.25 ;
disp_eeg(D1,offset,fs,ElecName, 'Data1' ) ;
saveas(gcf, 'Data1.png');

offset = max(max(abs(D2)))/2 ;
disp_eeg(D2,offset,fs,ElecName, 'Data2' ) ;
saveas(gcf, 'Data2.png');

N = 1024; %N_FFT
[F1,W1,K1] = COM2R(D1,21);
S1 = W1*D1;
for i = 1:21
    [psd1(i,:),w] = pwelch(S1(i,:),[],[],N,fs);
end
[F2,W2,K2] = COM2R(D2,21);
S2 = W2*D2;
for i = 1:21
    [psd2(i,:),w] = pwelch(S2(i,:),[],[],N,fs);
end

%Plotting Sources 
offset = max(max(abs(S1)))/2 ;
disp_eeg(S1,offset,fs,1:21, 'Data1: Sources' ) ;
ylabel('Sources','FontSize',10);
saveas(gcf, 'Data1;Source.png');

offset = max(max(abs(S2)))/2 ;
disp_eeg(S2,offset,fs,1:21, 'Data2: Sources' ) ;
ylabel('Sources','FontSize',10);
saveas(gcf, 'Data2;Source.png');

%Sources' PSD
offset = max(max(abs(psd1)))/2 ;
disp_eeg(psd1,offset,N/fs,1:21, 'Data1: Sources-PSD' ) ;
xlabel('Frequency (Hz)','FontSize',10);
ylabel('Sources','FontSize',10);
xlim([0,80]);
saveas(gcf, 'Data1;PSD.png');

offset = max(max(abs(psd2)))/2 ;
disp_eeg(psd2,offset,N/fs,1:21, 'Data2: Sources-PSD' ) ;
xlabel('Frequency (Hz)','FontSize',10);
ylabel('Sources','FontSize',10);
xlim([0,80]);
saveas(gcf, 'Data2;PSD.png');

%Sources' Spatial map
figure1 = figure('WindowState', 'maximized');
for i = 1:21
    subplot(3,7,i);
    plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F1(:,i));
    title(['Source ',num2str(i)]);
end
saveas(gcf, 'Data1;Space.png');


figure1 = figure('WindowState', 'maximized');
for i = 1:21
    subplot(3,7,i);
    plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F2(:,i));
    title(['Source ',num2str(i)]);
end
saveas(gcf, 'Data2;Space.png');

%Removing Noise Sources
Selsources1 = [3,5:9,11:21];
X_den1 = F1(:,Selsources1)*S1(Selsources1,:);
offset = max(max(abs(X_den1)))/1.2 ;
disp_eeg(X_den1,offset,fs,ElecName, 'Denoised Data1' ) ;
saveas(gcf, 'Denoised Data1.png');

Selsources2 = [4:6,8:13,15:21];
X_den2 = F2(:,Selsources2)*S2(Selsources2,:);
offset = max(max(abs(X_den2)))/1.2 ;
disp_eeg(X_den2,offset,fs,ElecName, 'Denoised Data2' ) ;
saveas(gcf, 'Denoised Data2.png');
