classdef Wall
    %class for walls on the map
    %walls must be horizontal or vertical only
    properties
        x1 %x coordinate of the start of the wall
        y1 %y coordinate of the start of the wall
        x2 %x coordinate of the end of the wall
        y2 %y coordinate of the end of the wall
        P1 %start wall point [x1;y1]
        P2 %end wall point [x2;y2], intersection with 'rue de la loi'
        n1 %normal vector positive of the wall
        n2 %normal vector negative of the wall
        street %street name
        eps_r = 5; %relative permittivity
    end
    methods
        function [obj] = Wall(x1,y1,x2,y2, street)
            obj.x1 = x1;
            obj.y1 = y1;
            obj.x2 = x2;
            obj.y2 = y2;
            obj.P1 = [x1;y1];
            obj.P2 = [x2;y2]; 
            obj.street = street;
            if(obj.x1==obj.x2)
                %vertical wall
                obj.n1=[1;0];
                obj.n2=[-1;0];
            end
            if(obj.y1==obj.y2)
                %horizontal wall
                obj.n1=[0;1];
                obj.n2=[0;-1];
            end
            plot(obj);
        end
        function plot(obj)
            plot([obj.x1 obj.x2],[obj.y1 obj.y2],'color','black','linewidth',2)
            hold on
        end
    end
end

