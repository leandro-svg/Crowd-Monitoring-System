function[antenna_image] = antenna_image(wall,point)
%return the position of the image antenna based on the wall and the point

    x1 = wall.Extremity1(1) ;
    y1 = round(wall.Extremity1(2)) ;
    x2 = wall.Extremity2(1);
    y2 = round(wall.Extremity2(2));
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