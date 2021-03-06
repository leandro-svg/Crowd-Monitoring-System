% warning('off', 'all')
% addpath functions data
% %% RT Simulator - 2022 MA2 Project
% 
% if exist('RT', 'var')
%     RT.delete();
% end
% 
% app = RTSimulatorMA2();

%% COMMUNICATION CHANNEL PART
%% Computation of the intersections
close all
 
setwindowLimits(app, [1,50,1,10]);
setTXLocations(app, 1, [2,5]); 
setRXLocations(app, 2, [8,3;23,5;38,9;48,2]);
RX = app.receivers2DPosition;
TX = app.emitters2DPosition;
paths = app.paths;
Walls = app.walls;
nbRealizations = 5;
density  = (0.01:0.01:3);
for j= 1:length(RX)
    filename1 = sprintf('RX_Position_%d.csv',RX(j,1)); 
    fclose(fopen(filename1, 'w'));
end

for dens =  (density)
    disp(dens)
    for p = 1:nbRealizations
        app.setPeopleDensity(dens);
        app.generatePeople(0);
        app.computeRays();
        people2DPosition = app.people2DPosition;
        P_wall = zeros(size(TX,1), size(RX,1), length(Walls),2);
        d = zeros(size(TX,1), size(RX,1), length(Walls),2);
        gamma_perp = zeros(size(TX,1), size(RX,1), length(Walls),1);
         for g = 1:size(TX,1)
            for j = 1:size(RX,1)
                for i = 1:length(Walls)
                    [antenna_image_TX] = antenna_image(Walls(i),TX(g,:)); 
                    %gamma_perp(i) = calculation_gamma_perp(Walls(i).nx, Walls(i).ny,antenna_image_TX, RX);
                    [P_wall(g, j, i,:), d(g, j, i,:)] = app.findIntersectionOnWall(antenna_image_TX', RX(j,:), Walls(i));
                    gamma_perp(g,j, i) = computeReflectionCoefficient_MIMO(app, Walls(i), d(g,j,i,:));
                end
            end
         end
        
        addpath("functions")
        %% Computation of the number of body crossed by the wall intersection ray
        nBC_wall = zeros(size(TX,1), size(RX,1),length(P_wall),1);
        for i = 1:size(TX,1)
            for j = 1:size(RX,1)
                for k=1:size(P_wall, 1)
                    nbScatterers = computeNumberOfBodiesCrossed_MIMO(TX(i, :), RX(j, :), P_wall(i,j,k,:),people2DPosition(1:end, :));
                    nBC_wall(i,j,k) = nbScatterers;
                end
            end
        end
        %% Computation of the CIR with the walls intersection
        %[A, R] = getCIR(rayInfo,1);
        [A_walls, R_walls] = getCIR_Walls_MIMO(TX,RX,rayInfo, P_wall,nBC_wall,0,gamma_perp);
        [Ntx, Nrx, Npeople] = size(rayInfo.paths);
        Npeople = Npeople-1;
        tau = R_walls./3e8;
        RLOS = R_walls - R_walls(:,1);
        tau = RLOS./3e8;
        min = -200;
        A_delimitited = min - 20*log10(abs(A_walls));
        
        B_list = [500e6];
        Bc = zeros(length(B_list), 1);
        tau_rms = zeros(length(B_list),1);
        mean_h = zeros(length(B_list),1);
        std_h = zeros(length(B_list),1);
        Path_loss =  zeros(length(B_list),1);
        
        for g = 1:size(TX,1)
            for j = 1:size(RX,1)
                for b = 1:length(B_list)
                    %% Creation of tapped delay line impulse response
                    B = B_list(b);
                    T = 1./B; %T is the delta tau, the width of a tap
                    %Creation of taps, which ray is going to which tap
                    indices = round((tau-tau(:,1))./T) + 1;%put every tap at zero / tau_LOS
                    Lmax = (max(indices(j,:)));
                    %number of taps
                    hFinite = zeros(Lmax, 1);
                    for i = 1:length(indices)
                        hFinite(indices(j,i)) = hFinite(indices(j,i)) + A_walls(j,i);%Put every ray in taps
                    end
                    delay_axis = (0:Lmax-1)*T;% Creation of the tap axis with delta tau times Lmax
                    
                    %
                    % RLOS = R - R(1);
                    %% Broadband attenuation
                    [Path_loss(b)] = broadband_attenuation(hFinite);
                    %disp(Path_loss(b))
                    filename = sprintf('RX_Position_%d.csv',RX(j,1)); 
                    writematrix([RX(j,:),B_list(b),dens, Npeople, Path_loss(b)], filename, 'WriteMode', 'append')
             
                    %% Delay spread
        %         
        %         %   [tau_rms(b), Bc(b)] = getDSandBC(hFinite);
        %             pdp = abs(hFinite.').^2; %power density
        %             Ptot = trapz(pdp); %total power
        %             tau_mean = trapz(delay_axis.*pdp) ./ Ptot; %mean tau
        %             tau_rms(b) = sqrt(trapz((delay_axis - tau_mean).^2.*pdp)./Ptot); %delay spread
        %             %%% problem with the tau inside the parenthesis
        %             Bc(b) = 1/(2*pi*tau_rms(b)); %whatever if we are in wide or narrow band
                end
            end
        end
    end
end
% Plot of path loss which respect to the number of people
close all
matrix_8 = readmatrix('RX_Position_8.csv');
matrix_23 = readmatrix('RX_Position_23.csv');
matrix_40 = readmatrix('RX_Position_38.csv');
matrix_58 = readmatrix('RX_Position_48.csv');
fclose(fopen('M2.csv', 'w'));


for i = 0:length(density)-1
    ans1  = mean(matrix_8(i*nbRealizations+1:(i+1)*nbRealizations,6));
    ans2  = mean(matrix_23(i*nbRealizations+1:(i+1)*nbRealizations,6));
    ans3  = mean(matrix_40(i*nbRealizations+1:(i+1)*nbRealizations,6));
    ans4  = mean(matrix_58(i*nbRealizations+1:(i+1)*nbRealizations,6));
    writematrix([matrix_58(i*nbRealizations+1,5),ans1,ans2, ans3, ans4], 'M2.csv', 'WriteMode', 'append')
end



% figure
% 
% stem(matrix(1:length(density),1), matrix(1:length(density),2))
% xlabel('Number of people');
% ylabel('Attenuation factor normalized');
% title('Correlation between the broadband attenuation and the number of people for receiver  at position 58'); 
% ylim([50,140])
% grid on
% 
% figure
% 
% stem(matrix(1:length(density),1), matrix(1:length(density),3))
% xlabel('Number of people');
% ylabel('Attenuation factor normalized');
% title('Correlation between the broadband attenuation and the number of people for receiver  at position 58'); 
% ylim([50,140])
% grid on
% 
% figure
% 
% stem(matrix(1:length(density),1), matrix(1:length(density),4))
% xlabel('Number of people');
% ylabel('Attenuation factor normalized');
% title('Correlation between the broadband attenuation and the number of people for receiver  at position 40'); 
% ylim([50,140])
% grid on
% figure
% 
% stem(matrix(1:length(density),1), matrix(1:length(density),5))
% xlabel('Number of people');
% ylabel('Attenuation factor normalized');
% title('Correlation between the broadband attenuation and the number of people for receiver  at position 8'); 
% ylim([50,140])
% grid on
density  = (0.01:0.01:3);
matrix = readmatrix('Neural_density_300samples.csv');
figure
e = std(matrix(1:length(density),2))*ones(size(matrix(1:length(density),1)));
errorbar(matrix(1:length(density),1),matrix(1:length(density),2),e)
title('position 8'); 
grid on

figure
e = std(matrix(1:length(density),3))*ones(size(matrix(1:length(density),1)));
errorbar(matrix(1:length(density),1),matrix(1:length(density),3),e)
title('position 23'); 
grid on


figure
e = std(matrix(1:length(density),4))*ones(size(matrix(1:length(density),1)));
errorbar(matrix(1:length(density),1),matrix(1:length(density),4),e)
title('position 40'); 
grid on


figure
e = std(matrix(1:length(density),5))*ones(size(matrix(1:length(density),1)));
errorbar(matrix(1:length(density),1),matrix(1:length(density),5),e)
title('position 58'); 
grid on

%% 

figure
X  = [matrix(:,1),matrix(:,2),matrix(:,3),matrix(:,4),matrix(:,5)];
Z = [matrix(:,1),matrix(:,2),matrix(:,3),matrix(:,4),matrix(:,5)];
maximum = max(matrix(:,1));
step_kNN = 0:round((maximum+1)/3):maximum+20;
Y = zeros(length(matrix(:,1)),1);
for i = 1:length(matrix(:,1))
    for j = 1:length(step_kNN)-1
        if (step_kNN(j) < matrix(i,1))
            if matrix(i,1) <= step_kNN(j+1) 
                Y(i) = step_kNN(j+1);
                X(i,1) = step_kNN(j+1);
            end
        end
    end
end
Mdl = fitcknn(X,Y);
[predictedLabels, posterior, cost]=predict(Mdl,Z);

trueLabels = Y;
size_posterior =  size(posterior);


confusionchart(trueLabels,predictedLabels,'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');






%%%%%%%%%%%%%%

