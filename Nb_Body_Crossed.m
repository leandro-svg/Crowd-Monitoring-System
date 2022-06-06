function [nBC_wall] = Nb_Body_Crossed(TXpos, RXpos, P_wall, path, people2DPosition)

nBC_wall = zeros(length(P_wall));
Npeople =  length(path)-1;
Ntx = size(TXpos, 1); % for all emitters;
Nrx = size(RXpos, 1);
app.paths = cell(Ntx, Nrx, Npeople);

%% Magic starts
for i = 1:Ntx
    for j = 1:Nrx
        for k=1:P_wall
            nbScatterers = app.computeNumberOfBodiesCrossed(TXpos(i, :), RXpos(j, :), P_wall(k, :),people2DPosition(1:end, :));
            % Compute intersections and distance
            intersectPt = [TXpos(i, :); P_wall(k, :); RXpos(j, :)];
            nBC_wall = nbScatterers;
        end
    end
end