% warning('off', 'all')
% addpath functions data
% %% RT Simulator - 2022 MA2 Project
% 
% if exist('RT', 'var')
%     RT.delete();
% end
% 
% app = RTSimulatorMA2();

%% COMMINCATION CHANNEL PART
% close all
% 
% setwindowLimits(app, [1,40,1,40]);
% setTXLocations(app, 1, [2,5]); 
% setRXLocations(app, 1, [38,2]);
% RX = app.receivers2DPosition;
% TX = app.emitters2DPosition;
% paths = app.paths;
% Walls = app.walls;
% nbRealizations = 10;
% density  = (0.01:0.3:6);
% [Ntx, Nrx, Npeople] = size(rayInfo.paths);
% Npeople = Npeople-1;
% error_A = zeros(length(density), nbRealizations);
% error_DS = zeros(length(density), nbRealizations);
% error_CB = zeros(length(density), nbRealizations);
% fclose(fopen('RX_app.csv', 'w'));
% firt=0;
% for dens =  (density)
%     firt = firt + 1;
%     disp(dens)
%     for p = 1:nbRealizations
%         app.setPeopleDensity(dens);
%         app.generatePeople(0);
%         app.computeRays();
%         people2DPosition = app.people2DPosition;
%         [A, R] = getCIR(rayInfo, 0);
%         
%         tau = R./3e8;
%         RLOS = R - R(1);
%         tau_CIR = RLOS./3e8;
%         min = -200;
%         A_delimitited = min - 20*log10(abs(A));
% %         figure;
% %         stem(tau_CIR*1e9, 20*log10(abs(A)), '*', 'LineWidth', 2, 'MarkerSize', 10, 'BaseValue',-200), grid on
% %         xlabel('tau [ns]');
% %         ylabel('|h(tau, t)| [dB]');
% %         title('Impulse Response (CIR)');
% %         
%         
%         B_list = [500e6];
%         Bc = zeros(length(B_list), 1);
%         tau_rms = zeros(length(B_list),1);
%         
% %         figure;
%         layout = tiledlayout(length(B_list), 2);
%         for b = 1:length(B_list)
%         
%             B = B_list(b);
%             T = 1./B; %T is the delta tau, the width of a tap
%         
%         
%             %Creation of taps, which ray is going to which tap
%             indices = round((tau-tau(1))./T) + 1;%put every tap at zero / tau_LOS
%             Lmax = max(indices);%number of taps
%             hFinite = zeros(Lmax, 1);
%             for i = 1:length(indices)
%                 hFinite(indices(i)) = hFinite(indices(i)) + A(i);%Put every ray in taps
%             end
%             delay_axis = (0:Lmax-1)*T;%Creation of the tap axis with delta tau times Lmax
%             %
%             % RLOS = R - R(1);
%             
%             [mean_h] = broadband_attenuation(hFinite);       
%             error_A(firt,p) = mean_h;
%             pdp = abs(hFinite.').^2; %power density
%             Ptot = trapz(pdp); %total power
%             tau_mean = trapz(delay_axis.*pdp) ./ Ptot; %mean tau
%             tau_rms(b) = sqrt(trapz((delay_axis - tau_mean).^2.*pdp)./Ptot); %delay spread
%             error_DS(firt, p) = tau_rms(b);
%             %%% problem with the tau inside the parenthesis
%             Bc(b) = 1/(2*pi*tau_rms(b)); %whatever if we are in wide or narrow band
%             error_CB(firt, p) = Bc(b);
%             
%         end
%     end
% end
% 
% % Plot of path loss which respect to the number of people



density  = (0.01:0.3:6);
density = density(1:13);
error_A =  error_A_full(1:13,:);
error_CB =  error_CB_full(1:13,:);
error_DS =  error_DS_full(1:13,:);
nbRealizations = 10;
figure
mean_PL = mean(error_A,2);
mean_PL = mean_PL';

std_vector = zeros(length(mean_PL),1);
for i = 1:length(error_A(:,1))
    std_vector(i) = std(error_A(i,:));
end
std_vector =  std_vector';
errorbar(density,mean_PL,std_vector)
title('Error bar of the path loss over crowd realization'); 
xlabel('Density [ppl/m^2]')
grid on

mean_DS = mean(error_DS,2);
mean_DS = mean_DS';
std_vector_DS = zeros(length(mean_PL),1);
for i = 1:length(error_DS(:,1))
    std_vector_DS(i) = std(error_DS(i,:));
end
std_vector_DS =  std_vector_DS';

mean_CB = mean(error_CB,2);
mean_CB = mean_CB';
std_vector_CB = zeros(length(mean_PL),1);
for i = 1:length(error_CB(:,1))
    std_vector_CB(i) = std(error_CB(i,:));
end
std_vector_CB =  std_vector_CB';


figure
subplot(3,2,1);
plot(density, mean_PL);
title('Mean path loss over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('mean of Path Loss [dBm]')
grid on

subplot(3,2,2);
plot(density, std_vector);
title('Standard deviation of the path loss over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('Standard deviation of Path Loss [dBm]')
grid on

subplot(3,2,3);
plot(density, mean_DS);
title('Mean of the delay spread over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('Mean of the delay spread [s]')
grid on

subplot(3,2,4);
plot(density, std_vector_DS);
title('Standard deviation of the delay spread over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('Std of the delay spread [s]')
grid on

subplot(3,2,5);
plot(density, mean_CB);
title('Mean of the coherence bandwidth over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('mean of Coherence bandwidth (Hz)')
grid on

subplot(3,2,6);
plot(density, std_vector_CB);
title('Standard deviation of the coherence bandwidth over crowd realization'); 
xlabel('Crowd Density [ppl/m^2]')
ylabel('std of Coherence bandwidth (Hz)')
grid on


%%%%%%%%%%%%%%

