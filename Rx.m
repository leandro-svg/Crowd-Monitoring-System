classdef Rx
    %class for vertical half-wave dipole receiver antenna object
    properties
        pos_x %coordinate in x of the position 
        pos_y %coordinate in y of the position 
        P %coordinate of the position as a point (x,y)      
    end
    methods
        function [obj] = Rx(pos_x,pos_y)
            obj.pos_x = pos_x;
            obj.pos_y = pos_y;
            obj.P = [pos_x;pos_y];
            %plot(obj);
        end
        function plot(obj)
            scatter(obj.pos_x,obj.pos_y,500,'magenta');
            hold on
        end
    end
end