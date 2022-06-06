function [A, R] = getCIR_Walls(TX, RX, rayInfo, P, nBC_wall, enableLOS, gamma_perp)
[Ntx, Nrx, Npath] = size(rayInfo.paths);
P_wall = P; % Position intersection on wall

BistaticConstant = 8.734510731023729e-08;
LOSConstant = 2.744027474528536e-06;
Beta = 51.3127;
Npath = Npath + length(P_wall(1,1,:,:));
h = zeros(Ntx, Nrx, Npath, 2);
h0 = zeros(Npath, 2);

for i = 1:Ntx
    for j = 1:Nrx
        for p = 1:Npath-size(P_wall(i,j,:,:),3)
            Rloss = prod(rayInfo.paths{i, j, p}.D.^2);
            Rphase = sum(rayInfo.paths{i,j,p}.D);
            h0(p,2) = Rphase;
            if ~rayInfo.paths{i,j,p}.type && enableLOS
                if rayInfo.paths{i,j,p}.nBC > 0
                    U = LOSConstant * 0.2524^rayInfo.paths{i,j,p}.nBC;
                else
                    U = LOSConstant;
                end
                h0(p,1) = sqrt(U/Rloss) * exp(1j*Beta*Rphase);
            elseif rayInfo.paths{i,j,p}.type == 1
                if rayInfo.paths{i,j,p}.nBC > 0
                    U = BistaticConstant * 0.2524^rayInfo.paths{i,j,p}.nBC;
                else
                    U = BistaticConstant;
                end
                h0(p, 1) = sqrt(U/Rloss) * exp(1j*Beta*Rphase);
            end

        end
        for p  = 1:1:size(P_wall(i,j,:,:),3)
            D_reflected = [P_wall(i,j,p,1),P_wall(i,j,p,2)];
            %D_reflected = [P_wall(p,1),P_wall(p,2)];
            r1 = sqrt(sum((TX(i, :) - D_reflected).^2)); % propagation distance from Tx to Human
            r2 = sqrt(sum((D_reflected - RX(j, :)).^2)); % propagation distance from Human to Rx
            
            D = [r1 r2]; % total distance traveled
            Rloss = prod(D.^2);
            Rphase = sum(D);

            U =  BistaticConstant * 0.2524^nBC_wall(i,j,p);
            h0(Npath-length(P_wall(i,j,:,:))+p,2) = Rphase;
            h0(Npath-length(P_wall(i,j,:,:))+p, 1) = gamma_perp(i,j,p)* sqrt(U/Rloss) * exp(1j*Beta*Rphase);
        end
        h(i, j, :, :) = sortrows(h0, 2);
    end
end
h = squeeze(h);
% [Ntx, Nrx, Npath, 2] = size(h);

% 
% R = real(h(:,:, 2));
% A = h(:,:, 1);

R = real(h(:,:, 2));
A = h(:,:,1);