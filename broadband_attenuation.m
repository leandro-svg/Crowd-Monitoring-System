function [Path_loss] = broadband_attenuation(hFinite)
P_transmitter = 0.01; %mW
E_tot = sum(hFinite);
P_receiver = abs(E_tot)^2;
Path_loss = P_transmitter/P_receiver;
Path_loss = 10*log10((Path_loss));
