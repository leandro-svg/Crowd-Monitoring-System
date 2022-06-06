%% DIGITAL COMMUNCATION PART

function [BER, SNR, MSE_freq, MSE_time] = digital_communication(h, plot_param,enableCIR,enableAWGN,enablePreamble,CIRknown, enableEqualization, enableFrequency, enableTime)

M = 16; %Modulation
cp_len =128; %Length of the cyclic prefix
SNR = (-5:1:30);
nbr_realization =20;
snr_realization = zeros(length(SNR),1);
mse_realization_freq = zeros(length(SNR),1);
mse_realization_time = zeros(length(SNR),1);

BER = zeros(length(SNR),1);
MSE_freq = zeros(length(SNR),1);
MSE_time = zeros(length(SNR),1);
h = h';
for k = 1:length(SNR)
    for c = 1:nbr_realization
%% Generation of the input data
        n_s= 20; %Number of symbol
       
        number_of_subcarriers = 1024; %Number of subcarrier channel
        Nb=n_s*number_of_subcarriers; %Total number of samples to be transmitted at the receiver
        Nbps = Nb/n_s; %Size of each OFDM block to add cyclic prefix

        
        % Generate random data source to be transmitted of length Nbps
        input_signal = randsrc(1, Nb, 0:M-1);
        if (enablePreamble)
            preamble = randsrc(1,2*number_of_subcarriers, 0:M-1);
            input_signal = [preamble,input_signal];
            n_s = n_s + 2;
            Nb=n_s*number_of_subcarriers;
            Nbps = Nb/n_s;
        end
        if (plot_param)
            figure;stem(input_signal); grid on; xlabel('Data Points'); ylabel('Amplitude')
            title('Original Data ')
        end
        % Perform QPSK modulation on the input source data
        D_tx = (qammod(input_signal, M));
        if (plot_param)
            figure;stem(D_tx);title('QPSK Modulation ')
        end
        % Generate the data matrix
        D_tx_matrix = reshape(D_tx, number_of_subcarriers, n_s);

...........................................................................
%% IFFT AND CYCLIC PREFIX
        ifft_D_tx = ifft(D_tx_matrix, number_of_subcarriers,1);
         S_tx_prefix_Added = [ifft_D_tx(end-cp_len+1:end,:);ifft_D_tx];

        
        
        [w_stx, l_stx]=size(S_tx_prefix_Added);
        len_stx = w_stx*l_stx;
        S_tx_cp = reshape(S_tx_prefix_Added, 1, len_stx);
        if (plot_param)
            figure,plot(real(S_tx_cp)); xlabel('Time'); ylabel('Amplitude');
            title('OFDM Signal');grid on;
        end
...........................................................................        
 %% CIR
        alpha = rand(1)*0.5;
        beta = sqrt(0.5^2-alpha^2);
        
        h =  [1 0 alpha  + 1j*beta];
        if (enableCIR)
            S_channel = conv(S_tx_cp,h);
            S_channel = S_channel(1:numel(S_tx_cp));
        else 
            S_channel = S_tx_cp;
        end
 %% CONVOLVE NOISE
        if (enableAWGN)
            S_rx_cp = awgn(S_channel, SNR(k), 'measured');
        else
            S_rx_cp = S_channel;
        end
...........................................................................
 
 %% FOURIER TRANSFORM
        
        S_rx_matrix = reshape(S_rx_cp,  (number_of_subcarriers+cp_len), n_s);
        D_rx_fft = fft(S_rx_matrix(cp_len+1:end,:) , number_of_subcarriers, 1);
        OFDM_signal = reshape(D_rx_fft, 1,Nb);
 %% CHANNEL TRANSFER FONCTION FREQUENCY DOMAIN
        if (CIRknown)
              H_est_freq = fft(h, number_of_subcarriers);
        else
              preamble_x_1 = D_tx(1:number_of_subcarriers);   % Pilot extraction
              preamble_y_1 = OFDM_signal(1:number_of_subcarriers);   % Pilot extraction from the DQPSK received data
              H_est_1 = (preamble_y_1./preamble_x_1);
              preamble_x_2 = D_tx(number_of_subcarriers+1:2*number_of_subcarriers);   % Pilot extraction
              preamble_y_2 = OFDM_signal(number_of_subcarriers+1:2*number_of_subcarriers);   % Pilot extraction from the DQPSK received data
              H_est_2 = (preamble_y_2./preamble_x_2);
              H_est_freq = (H_est_1+H_est_2)/2;
              
        end
        H_CIR = fft(h, number_of_subcarriers);
        if (plot_param)
            figure
            H_CIR = fft(h, number_of_subcarriers);
            plot(1:number_of_subcarriers,abs(H_CIR))
            hold on
            plot(1:number_of_subcarriers, abs(H_est_freq))
            title('TF')
        end


 ...........................................................................
%% CHANNEL TRANSFER FONCTION TIME DOMAIN OR FREQUENCY DOMAIN
        
        h_est = ifft(H_est_freq);
        %semilogy(1:number_of_subcarriers, abs(h_est).^2);
        H_est_time = fft(h_est(1:cp_len), number_of_subcarriers);

        if(enableTime)
            H_est = H_est_time;
        elseif(not(enableTime))
            H_est = H_est_freq;
        end     
 ...........................................................................       
%% COMPUTATION MSE BETWEEN H KNOWN IMPULSE RESPONSE AND H ESTIMATED
        mse_H_frequency = MSE(H_CIR, H_est_freq);
        mse_H_time = MSE(H_CIR, H_est_time);
        mse_realization_freq(k) = mse_H_frequency + mse_realization_freq(k);
        mse_realization_time(k) = mse_H_time + mse_realization_time(k);
 ...........................................................................
%% EQUALIZATION
        OFDM_equalizer = zeros(size(OFDM_signal));
        if (enableEqualization)
            for v = 1:n_s
               OFDM_equalizer((v-1)*number_of_subcarriers+1:v*number_of_subcarriers) = OFDM_signal((v-1)*number_of_subcarriers+1:v*number_of_subcarriers)./H_est; %ZF equalizer
            end
        else
            OFDM_equalizer = OFDM_signal;
        end
...........................................................................
%% DEMODULATION
        output_signal = qamdemod(OFDM_equalizer,M);
        if (enablePreamble)
            output_signal = output_signal(2:end); %Remove the preambles
            input_signal = input_signal(2:end);
        end
        if (plot_param)
            figure
            stem(input_signal)
            hold on
            stem(output_signal,'rx');
            grid on;xlabel('Data Points');ylabel('Amplitude');
            title('Recieved Signal with error')
        end
        [~,ratio] = biterr(de2bi(input_signal,(M)), de2bi(output_signal, (M)));
        snr_realization(k) = snr_realization(k) + ratio;
        
 ...........................................................................
    end
    BER(k) = snr_realization(k)/nbr_realization;
    MSE_freq(k) = mse_realization_freq(k)/nbr_realization;
    MSE_time(k) = mse_realization_time(k)/nbr_realization;
end





