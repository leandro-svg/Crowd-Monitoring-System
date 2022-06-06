function app = getRays(app, TXpos, RXpos)

wb = waitbar(0, 'Computing the rays...');


progression = 0;
Ntx = size(TXpos, 1); % for all emitters;
Nrx = size(RXpos, 1);
N_calculations = Npeople*Ntx*Nrx; % number of calculations, i.e, number_of_people * number_of_TX * number_of_RX

app.paths = cell(Ntx, Nrx, Npeople);

%% Magic starts
for i = 1:Ntx
    for j = 1:Nrx
        ind = 1;
        %% Direct path computation
        if app.param.AllowDirectPath.Value % if direct path is allowed from GUI
            [~, nbScatterers] = computeNormalVectorAndNbScatterers(app, TXpos(i, :), RXpos(j, :), app.people2DPosition);
            
            app.paths{i,j,ind}.type = 0; % type
            app.paths{i,j,ind}.intersectPt = [TXpos(i, :); RXpos(j, :)];
            app.paths{i,j,ind}.D = sqrt(sum((TXpos(i, :) - RXpos(j, :)).^2)); % total propagation distance of the direct path          
            temp_ang = RXpos(j, :) - TXpos(i, :);
            app.paths{i,j,ind}.angularInfo = atand(temp_ang(2)/temp_ang(1));
            app.paths{i,j,ind}.nBC = nbScatterers;
            ind = ind + 1;
        end
        
        %% Reflected path computation (only single bounces)
        for k=1:Npeople
            progression = progression + 1;
            waitbar(progression/N_calculations, wb, 'Computing the rays...');
            
            nbScatterers = computeNumberOfBodiesCrossed(TXpos(i, :), RXpos(j, :), app.people2DPosition(k, :), app.people2DPosition(1:end ~= k, :));
            % Compute intersections and distance
            app.paths{i,j,ind}.type = 1;
            app.paths{i,j,ind}.intersectPt = [TXpos(i, :); app.people2DPosition(k, :); RXpos(j, :)];
            
            r1 = sqrt(sum((TXpos(i, :) - app.people2DPosition(k, :)).^2)); % propagation distance from Tx to Human
            r2 = sqrt(sum((app.people2DPosition(k, :) - RXpos(j, :)).^2)); % propagation distance from Human to Rx
            
            app.paths{i,j,ind}.D = [r1 r2]; % total distance traveled
            temp_ang1 = app.people2DPosition(k, :) - TXpos(i,:);
            temp_ang2 = RXpos(j,:) - app.people2DPosition(k, :);
            app.paths{i,j,ind}.angularInfo = [atand(temp_ang1(2)/temp_ang(1)) atand(temp_ang2(2)/temp_ang2(1))];
            app.paths{i,j,ind}.nBC = nbScatterers;
            ind = ind + 1;
        end
    end
end
waitbar(1, wb, 'Finished');
close(wb);

end



%% Computations below are just geometric stuff to obtain some terms, can be ignored
function nbScatterers = computeNumberOfBodiesCrossed(TXpos, RXpos, people2DPositionA, people2DPositionB)

% Incident ray
nbScatterers1 = computeBC(TXpos, people2DPositionA, people2DPositionB);

%Reflected ray
nbScatterers2 = computeBC(people2DPositionA, RXpos, people2DPositionB);

nbScatterers = nbScatterers1 + nbScatterers2;

end

function nBC = computeBC(startPoint, endPoint, peoplePoints)
nBC = 0;
distance_finder = @(point1, point2) sqrt(sum((point1 - point2).^2));
lineDistance = distance_finder(startPoint, endPoint);
for i = 1:size(peoplePoints, 1)
    currentPerson = peoplePoints(i, :);
    start2point = distance_finder(startPoint, currentPerson);
    point2end = distance_finder(currentPerson, endPoint);
    if abs(lineDistance - (start2point + point2end)) < 0.01
        nBC = nBC + 1;
    end
end
% nBC = 1;
end
function [n, nbScatterers] = computeNormalVectorAndNbScatterers(app, startPoint, endPoint, people)
% Compute the unitary vector n between startPoint and endPoint.
% Then, for each person in people, look if he/she is not
% crossed. The number of people crossed in given in
% nbScatterers.

P12 = endPoint - startPoint; % Vector from startPoint to endPoint
n = P12 / norm(P12);  % Normalized vector in direction of the ray

% Seach for intersection with people
if ~isempty(people)
    startToPerson = people - startPoint; % Line from startPoint to the people center
    v = abs(n(1) * startToPerson(:,2) - n(2) * startToPerson(:,1)); % Norm of the 2D cross-product
    doIntersect = (v <= app.peopleRadius); % Shortest distance between ray and center of the disc
    nbScatterers = sum(doIntersect);
else
    nbScatterers = 0;
end
end
