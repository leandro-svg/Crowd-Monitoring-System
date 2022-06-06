function [mean_h, std_h] = attenuation_factor(h)
h = h';
mean_h = abs(mean(h));
std_h = std(h);