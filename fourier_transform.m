function [H,t_axis] = fourier_transform(h_0, B)
Q = 1024;
t_axis = [-Q/2:Q/2]*(B/Q);
H = fft(h_0, Q);
