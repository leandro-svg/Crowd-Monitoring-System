warning('off', 'all')
addpath functions data
%% RT Simulator - 2022 MA2 Project

if exist('RT', 'var')
    RT.delete();
end

RTSimulatorMA2();

%% COMMINCATION CHANNEL PART
close all
[A, R] = getCIR(rayInfo, 0);
tau = R./3e8;
RLOS = R - R(1);
tau_CIR = RLOS./3e8;
min = -200;
A_delimitited = min - 20*log10(abs(A));
figure;
stem(tau_CIR*1e9, 20*log10(abs(A)), '*', 'LineWidth', 2, 'MarkerSize', 10, 'BaseValue',-200), grid on
xlabel('tau [ns]');
ylabel('|h(tau, t)| [dB]');
title('Impulse Response (CIR)');





B_list = [10e6, 50e6, 100e6, 500e6, 1e9];
Bc = zeros(length(B_list), 1);
attenuation_factor = zeros(length(B_list), 1);
tau_rms = zeros(length(B_list),1);

figure;
layout = tiledlayout(length(B_list), 2);
for b = 1:length(B_list)

    B = B_list(b);
    T = 1./B; %T is the delta tau, the width of a tap

    [H, t_axis] = fourier_transform(A, B);


    %Creation of taps, which ray is going to which tap
    indices = round((tau-tau(1))./T) + 1;%put every tap at zero / tau_LOS
    Lmax = max(indices);%number of taps
    hFinite = zeros(Lmax, 1);
    for i = 1:length(indices)
        hFinite(indices(i)) = hFinite(indices(i)) + A(i);%Put every ray in taps
    end
    delay_axis = (0:Lmax-1)*T;%Creation of the tap axis with delta tau times Lmax
    %
    % RLOS = R - R(1);
    
    [Path_loss(b)] = broadband_attenuation(hFinite);
    pdp = abs(hFinite.').^2; %power density
    Ptot = trapz(pdp); %total power
    tau_mean = trapz(delay_axis.*pdp) ./ Ptot; %mean tau
    tau_rms(b) = sqrt(trapz((delay_axis - tau_mean).^2.*pdp)./Ptot); %delay spread
    %%% problem with the tau inside the parenthesis
    Bc(b) = 1/(2*pi*tau_rms(b)); %whatever if we are in wide or narrow band
    
    h=nexttile;
    [H, t_axis] = fourier_transform((hFinite), B);
    plot(t_axis(1:end-1), 20*log10(abs(H))), grid on, hold on
    xlabel('Frequency [MHz]');
    ylabel('|H(f)| [dB]');
    title('Transfer function'); 
    
    h=nexttile;

    stem(delay_axis*1e9, 20*log10(abs(hFinite)), '*', 'LineWidth', 2, 'MarkerSize', 10, 'BaseValue',-200), grid on, hold on
    lim =  max(delay_axis);
    %xlim([0 40])
%     keyboard    
    %h.XAxis.TickLabels = {'0','20','40', '60', '80', '100'};
    ylabel('PDP [dB]'), xlabel('Delay [ns]')
    title(['BW = ' num2str(B/1e6) 'MHz Delay spread = ' num2str(tau_rms(b)*1e9) 'ns ' 'Coh. Bw. = ' num2str(Bc(b)*1e-6) 'MHz'], 'Interpreter','latex')

%     title(layout,'Fading variability')
%     h= nexttile;
%     faxis = (-512:511)*B/1024;
%     plot(faxis*1e-6, 20*log10(abs(fft(hFinite, 1024)))), grid on
%     xlabel('Frequency [MHz]'), ylabel('abs(H) [dB]')
%     title(['BW = ' num2str(B/1e6) 'MHz Delay spread = ' num2str(tau_rms(b)*1e9) 'ns ' 'Coh. Bw. = ' num2str(Bc(b)*1e-6) 'MHz'], 'Interpreter','latex')
%     %title(['BW = ' num2str(B/1e6) 'MHz Delay spread = ' num2str(delay_spread*1e9) 'ns'], 'Interpreter','latex')
%     stem(RLOS/3e8, abs(A))
    %
end

%%%%%%%%%%%%%%

