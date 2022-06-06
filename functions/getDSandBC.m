% function [tau_rms, Bc] = getDSandBC(finiteCIR, delay_axis)
%% this is sort of a pseudo code
% finiteCIR: a vector containing the complex amplitudes per tap
% delay_axis: a vector containing the tau values per entry of finiteCIR
% make sure that finiteCIR is sorted wrt. delay_axis

pdp = abs(finiteCIR).^2;
P = trapz(pdp);
tau_m = 1/P*trapz(delay_axis.*pdp);
tau_rms = sqrt(1/P*trapz((delay_axis).^2 .* pdp) - tau_m.^2);
Bc = 1/(2*pi*tau_rms);

% end