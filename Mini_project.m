warning('off', 'all')
addpath functions data
%% RT Simulator - 2022 MA2 Project

if exist('RT', 'var')
    RT.delete();
end

app = RTSimulatorMA2();

%% COMMUNICATION CHANNEL PART
%% Computation of the intersections
close all
fclose(fopen('M.csv', 'w'));
setwindowLimits(app, [1,60,1,10]);
setTXLocations(app, 1, [2,5]);
setRXLocations(app, 1, [58,4]);
RX = app.receivers2DPosition;
TX = app.emitters2DPosition;
paths = app.paths;
Walls = app.walls;
nbRealizations = 1;
density  = [0.01:0.01:1];
for dens =  (density)
    disp(dens)
    for p = 1:nbRealizations
        app.setPeopleDensity(dens);
        app.generatePeople(0);
        app.computeRays();
        people2DPosition = app.people2DPosition;
        P_wall = zeros(length(Walls),2);
        nBC_wall = zeros(length(P_wall),1);
        d = zeros(length(Walls),2);
        gamma_perp = zeros(length(Walls),1);
        for i = 1:length(Walls)
            [antenna_image_TX] = antenna_image(Walls(i),TX); 
            %gamma_perp(i) = calculation_gamma_perp(Walls(i).nx, Walls(i).ny,antenna_image_TX, RX);
            [P_wall(i,:), d(i,:)] = app.findIntersectionOnWall(antenna_image_TX', RX, Walls(i));
            gamma_perp(i) = computeReflectionCoefficient(app, Walls(i), d(i,:));
        end
        
        addpath("functions")
        %% Computation of the number of body crossed by the wall intersection ray
        for i = 1:size(TX,1)
            for j = 1:size(RX,1)
                for k=1:size(P_wall, 1)
                    nbScatterers = computeNumberOfBodiesCrossed(TX(i, :), RX(j, :), P_wall(k, :),people2DPosition(1:end, :));
                    nBC_wall(k) = nbScatterers;
                end
            end
        end
        %% Computation of the CIR with the walls intersection
        %[A, R] = getCIR(rayInfo,1);
        [A_walls, R_walls] = getCIR_Walls(TX, RX,rayInfo, P_wall,nBC_wall,0, gamma_perp);
        [Ntx, Nrx, Npeople] = size(rayInfo.paths);
        Npeople = Npeople-1;
        tau = R_walls./3e8;
        RLOS = R_walls - R_walls(1);
        tau = RLOS./3e8;
        min = -200;
        A_delimitited = min - 20*log10(abs(A_walls));
    %     figure;
    %     stem(tau*1e9, 20*log10(abs(A_walls)), '*', 'LineWidth', 2, 'MarkerSize', 10, 'BaseValue',-200), grid on
    %     xlabel('tau [ns]');
    %     ylabel('|h(tau, t)| [dB]');
    %     title('Impulse Response (CIR)');
        
        B_list = [500e6];
        Bc = zeros(length(B_list), 1);
        tau_rms = zeros(length(B_list),1);
        mean_h = zeros(length(B_list),1);
        std_h = zeros(length(B_list),1);
        Path_loss =  zeros(length(B_list),1);
        
        %figure;
        layout = tiledlayout(length(B_list),1);
        for b = 1:length(B_list)
            %% Fourier transorm of IR
            B = B_list(b);
            T = 1./B; %T is the delta tau, the width of a tap
        
            [H, t_axis] = fourier_transform(A_walls, B);
        
            %% Creation of tapped delay line impulse response
        
            %Creation of taps, which ray is going to which tap
            indices = round((tau-tau(1))./T) + 1;%put every tap at zero / tau_LOS
            Lmax = max(indices);%number of taps
            hFinite = zeros(Lmax, 1);
            for i = 1:length(indices)
                hFinite(indices(i)) = hFinite(indices(i)) + A_walls(i);%Put every ray in taps
            end
            delay_axis = (0:Lmax-1)*T;% Creation of the tap axis with delta tau times Lmax
            %
            % RLOS = R - R(1);
            %% Broadband attenuation
            [Path_loss(b)] = broadband_attenuation(hFinite);
            %disp(Path_loss(b))
            
        
            %% Delay spread
        
        %   [tau_rms(b), Bc(b)] = getDSandBC(hFinite);
            pdp = abs(hFinite.').^2; %power density
            Ptot = trapz(pdp); %total power
            tau_mean = trapz(delay_axis.*pdp) ./ Ptot; %mean tau
            tau_rms(b) = sqrt(trapz((delay_axis - tau_mean).^2.*pdp)./Ptot); %delay spread
            %%% problem with the tau inside the parenthesis
            Bc(b) = 1/(2*pi*tau_rms(b)); %whatever if we are in wide or narrow band
        
            writematrix([B_list(b),dens, Npeople, Path_loss(b), tau_rms(b)], 'M.csv', 'WriteMode', 'append')
            %% Plot of Transfer function and TDL IR
        
    %         h=nexttile;
    %         semilogy(B_list, mean_h(b)), grid on, hold on
    %         xlabel('Frequency [MHz]');
    %         ylabel('Attenuation factor');
    %         title('Attenuation Factor'); 
        %     title(layout,'Fading variability')
        
        
        
            
        end
    end
end
% Plot of path loss which respect to the number of people
close all
fclose(fopen('M2.csv', 'w'));
matrix = readmatrix('M.csv');

for i = 0:length(density)-1
    ans  = mean(matrix(i*nbRealizations+1:(i+1)*nbRealizations,4));
    ans2 = mean(matrix(i*nbRealizations+1:(i+1)*nbRealizations,5));
    writematrix([matrix(i*nbRealizations+1,3),ans, ans2], 'M2.csv', 'WriteMode', 'append')
end

matrix = readmatrix('M2.csv');
figure

stem(matrix(:,1), matrix(:,2))
xlabel('Number of people');
ylabel('Attenuation factor normalized');
title('Correlation between the broadband attenuation and the number of people'); 
ylim([60,140])
grid on


figure
e = std(matrix(:,2))*ones(size(matrix(:,1)));
errorbar(matrix(:,1),matrix(:,2),e)
grid on

%% 

figure
X  = [matrix(:,1),matrix(:,2),matrix(:,3)];
Z = [matrix(:,1),matrix(:,2),matrix(:,3)];
maximum = max(matrix(:,1));
step_kNN = 0:round((maximum+1)/5):maximum+20;
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

