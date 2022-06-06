
% warning('off', 'all')
% addpath functions data
% %% RT Simulator - 2022 MA2 Project
% 
% if exist('RT', 'var')
%     RT.delete();
% end
% 
% RTSimulatorMA2();
%% DIGITALCOMMINCATION PART
%close all

[h, R] = getCIR(rayInfo, 0);
plot_param = 0;
enableTime = 1;
figure
%% OFDM without noise and CIR
% disp('OFDM without noise and CIR')
% enableCIR = 0;
% enableAWGN = 0;
% enablePreamble = 0;
% CIRknown = 1;
% enableEqualization = 0;
% [BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization,enableFrequency, enableTime);
% semilogy(SNR, BER); 
% 
% hold on
% %% OFDM with noise without CIR
% disp('OFDM noise')
% enableCIR = 0;
% enableAWGN = 1;
% enablePreamble = 0;
% CIRknown = 1;
% enableEqualization = 0;
% [BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization,enableFrequency, enableTime);
% semilogy(SNR, BER); 
% hold on

% %% OFDM with noise without equalization
% disp('OFDM no equalization')
% enableCIR = 1;
% enableAWGN = 1;
% enablePreamble = 0;
% CIRknown = 1;
% enableEqualization = 0;
% [BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization, enableFrequency, enableTime);
% semilogy(SNR, BER); 
% hold on



%% OFDM with noise and known CIR with equalization
% disp('OFDM known CIR with equalization')
% enableCIR = 1;
% enableAWGN = 1;
% enablePreamble = 0;
% CIRknown = 1;
% enableEqualization = 1;
% [BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization, enableFrequency, enableTime);
% semilogy(SNR, BER, '--x'); 
% hold on

%% OFDM with noise and unknown CIR, preamble added with equalization in time domain
disp('OFDM unknown CIR with preamble with equalization in time domain')
enableCIR = 1;
enableAWGN = 1;
enablePreamble = 1;
CIRknown = 0;
enableEqualization = 1;
[BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization,enableFrequency, enableTime);
semilogy(SNR, BER, '--o'); 
% hold on

% %% OFDM with noise and unknown CIR, preamble added with equalization in frequency domain
% disp('OFDM unknown CIR with preamble with equalization in frequency domain')
% enableCIR = 1;
% enableAWGN = 1;
% enablePreamble = 1;
% CIRknown = 0;
% enableEqualization = 1;
% enableTime = 0;
% [BER, SNR, MSE_freq, MSE_time] = digital_communication(h,plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization,enableFrequency, enableTime);
% semilogy(SNR, BER, '--x'); 
legend('OFDM unknown CIR with preamble with equalization in time domain')
%legend('OFDM noisy','OFDM no equalization','OFDM known CIR with equalization',  'OFDM unknown CIR with preamble with equalization in time domain','OFDM unknown CIR with preamble with equalization in freqeuncy domain')
grid on;xlabel('SNR (dB)'); ylabel('BER')
title('OFDM BER ')

figure
plot(SNR, MSE_freq, '--o');
hold on
plot(SNR, MSE_time, '--o'); 
legend('H et estimated H in Frequency domain channel estimation','H and estimated H in Time domain channel estimation')
grid on;xlabel('SNR (dB)'); ylabel('MSE')
title('MSE')


% 
% SNR = (0:1:20);
% Pb = zeros(length(SNR),1);
% for  i = 1:length(SNR)
%     Pb(i) = 0.5*(1-sqrt(SNR(i)/(SNR(i)+1)));
% end
% semilogy(SNR, Pb); 

