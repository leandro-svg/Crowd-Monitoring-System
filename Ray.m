classdef Ray
    %class for rays that goes from tx to rx
    properties
        gamma = 1 %coef reflection 
        D = 1 %coef diffraction
        d %distance traveled
        LOS = 0 % ==1 if the ray is from LOS
        tau %propagation delay
        ground_R = 0 % == 1 if the ray has ground reflection
        theta_i_ground %incidence angle with the ground reflection
        nbr_R = 0 %number of reflection 
        theta_i = 0 % sum of incidence angle 
    end
    methods
        function [obj] = Ray(antenna_image,point_rx)  
            c = 3*10^8; %speed of light
            obj.d = norm(point_rx - antenna_image);
            obj.tau = obj.d/c;
        end
    end
end