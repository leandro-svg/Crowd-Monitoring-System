function[rays] = find_rays(wall_list,tx,rx)
%return a list with all the rays coming to rx from tx considering all the walls 
%in the map
%draw = 1 to draw the rays
    
    rays = {}; %cell array containing all the objects ray
    %% LOS
    [ray_LOS,LOS] = line_of_sight(wall_list,tx,rx); %check if there is a LOS ray
    if LOS
        ray_LOS.LOS = 1;
        rays{end+1} = ray_LOS;
    end

    %% Ground reflection  
    if LOS
        [ray_ground] = ground_reflexion(rx.P,tx.P,draw);       
        rays{end+1} = ray_ground;
    end
    
    %% Wall reflection for rx in "rue de la loi"
    
        %% One reflection
    for k=1:length(wall_list)
        wall = wall_list(k);
        antenna_image1 = antenna_image(wall,tx.P);% calculate the image antenna of the wall wrt to the position of the transmitter 
        [impact1,pt_impact1] = Impact(wall,antenna_image1,rx.P);% check if the line segment from antenna image and rx intersects the wall
        if impact1 % if true, there might be a ray with 1 reflection if that ray doesnt intersect another wall
            intersection_1a = check_intersection_wall_list(wall_list,tx.P,pt_impact1); % check if the line segment from tx to the reflecting point intersects another wall
            intersection_1b = check_intersection_wall_list(wall_list,pt_impact1,rx.P); % check if the line segment from the reflecting point to rx intersects another wall
            intersection_1 = intersection_1a + intersection_1b; %% == 0 if no intersection for all line segment of the ray
            if not(intersection_1) %there is a ray with 1 reflection
                ray = Ray(antenna_image1,rx.P);%create a ray
                [theta_i] = calculation_theta_i(wall,antenna_image1,rx.P);%calculate the incidence angle on the reflecting wall
                gamma_perp = calculation_gamma_perp(theta_i);%calculate the reflection coef
                ray.gamma = ray.gamma * gamma_perp;%change the reflection coef of that ray
                ray.theta_i = ray.theta_i + theta_i;
                ray.nbr_R = 1;  
                if draw
                    plot([tx.pos_x pt_impact1(1) rx.pos_x],[tx.pos_y pt_impact1(2) rx.pos_y])%draw the ray
                end      
                rays{end+1} = ray;
            end
        end
    end
end

function[theta_i] = calculation_theta_i(wall,point1,point2)
%Calculate the incidence angle of the line composed of point1 and point2 on
%the wall
    v = point2 - point1;
    ang1 = acosd(dot(wall.n1,v)/(norm(wall.n1)*norm(v)));%angle entre rayon incident et normale positive
    ang2 = acosd(dot(wall.n2,v)/(norm(wall.n2)*norm(v)));%angle entre rayon incident et normale negative
    theta_i = min(ang1,ang2);
end


function[ray_LOS,LOS] = line_of_sight(wall_list,tx,rx,draw)
    %check is there is a LOS between tx and rx, if yes: create a ray
    ray_LOS = [];
    LOS = 1;
    for wall = wall_list
        [impact,pt_impact] = Impact(wall,tx.P,rx.P);
        if impact 
        %If there is at least 1 intersection between a wall and the ray
        %between rx and tx, there is no LOS
            LOS = 0; 
        end
    end
    if LOS
        ray_LOS = Ray(tx.P,rx.P);%create a ray
        if draw
            plot([tx.P(1) rx.P(1)],[tx.P(2) rx.P(2)]);
        end
    end
end     

function [intersection] = check_intersection_wall_list(wall_list,point1,point2)
% Check if there is at least 1 intersection between at least one wall from wall_list
% and the ray between point1 and point2
    intersection = 0;
    for wall = wall_list
                [impact,pt_impact] = Impact(wall,point1,point2);
        if impact 
            intersection = 1; 
        end
    end
end

function[impact,pt_impact] = Impact(wall,point1,point2)
%check if there is an intersection between the line segment of point1 and 
%point 2 and the wall
%if yes: impact == 1 and pt_impact == point of intersection

    x1 = wall.x1 ;
    y1 = wall.y1 ;
    x2 = wall.x2;
    y2 = wall.y2;

    x3 = point1(1);
    y3 = point1(2);
    v1 = (point2(1) - x3);
    v2 = (point2(2) - y3);   
    x4 = x3 + v1;
    y4 = y3 + v2;
    pt_impact = [NaN;NaN];
    
    den = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
    if (den == 0)
        impact = 0;
    else
        t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/den;
        u = -((x1-x2)*(y1-y3)-(y1-y2)*(x1-x3))/den;
        if(t>0 && t<1 && u>0 && u<1 )
            impact = 1;
            impact_x = x1 + t*(x2-x1);
            impact_y = y1 + t*(y2-y1);
            pt_impact = [impact_x;impact_y];
            %plot([impact_x rx_x],[impact_y rx_y]);
            %scatter(impact_x,impact_y);
        else
            impact = 0;
        end
    end

end

function[antenna_image] = antenna_image(wall,point)
%return the position of the image antenna based on the wall and the point

    x1 = wall.x1 ;
    y1 = wall.y1 ;
    x2 = wall.x2;
    y2 = wall.y2;
    if(x1==x2)        % if the wall is vertical 
         antenna_image(1) = 2*x1-point(1); 
         antenna_image(2) = point(2) ;               
    elseif(y1 == y2)                      % if the wall is horizontal
         antenna_image(1) = point(1);  
         antenna_image(2) = 2*y1-point(2);
    end
    
    %scatter(antenne_image(1),antenne_image(2));
    antenna_image = [antenna_image(1);antenna_image(2)];
end

function [ray_ground] = ground_reflexion(rx_pos,tx_pos,draw)
    global c
    ray_ground = Ray(tx_pos,rx_pos);
    h = 2; %height of the user equipements (m)
    
    LOS_d = norm(rx_pos - tx_pos); %line of sight distance
    distance = 2*sqrt(h^2+(LOS_d/2)^2);
    ray_ground.d = distance;
    ray_ground.tau = distance/c;
    
    theta_i = atand(LOS_d/(2*h)); %degree
    gamma_para = calculation_gamma_para(theta_i);
    ray_ground.gamma = ray_ground.gamma * gamma_para;
    ray_ground.ground_R = 1;
    ray_ground.theta_i_ground = theta_i;
    ray_ground.theta_i = theta_i;
    pt_reflection = [(rx_pos(1)+tx_pos(1))/2;(rx_pos(2)+tx_pos(2))/2];
    
    if(draw)
        scatter(pt_reflection(1),pt_reflection(2),300,'x');
        hold on
    end
end


function [gamma_para] = calculation_gamma_para(theta_i)
%calculation of reflection coefficient for parallel polarization based on the relative permittivity 
%and the angle of incidence
global epsilon_r
num = cosd(theta_i) - (1/sqrt(epsilon_r))*sqrt(1-(sind(theta_i)^2/epsilon_r));
den = cosd(theta_i) + (1/sqrt(epsilon_r))*sqrt(1-(sind(theta_i)^2/epsilon_r));
gamma_para = num/den;
end

function [gamma_perp] = calculation_gamma_perp(theta_i)
%calculation of reflection coefficient for perpendicular polarization based on the relative permittivity 
%and the angle of incidence
global epsilon_r
num = cosd(theta_i) - sqrt(epsilon_r)*sqrt(1-(sind(theta_i)^2/epsilon_r));
den = cosd(theta_i) + sqrt(epsilon_r)*sqrt(1-(sind(theta_i)^2/epsilon_r));
gamma_perp = num/den;
end


