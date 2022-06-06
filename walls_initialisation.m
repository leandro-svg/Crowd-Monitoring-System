function[wall_list] = walls_initialisation()
%create a list containing all the objects of class Wall on the map



    m1 = Wall(0,10,10,10, "loi");
    m2 = Wall(0,0,10,0, "loi");
    
  
    wall_list = [m1,m2];
end

