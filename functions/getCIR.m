function [A, R] = getCIR(rayInfo, enableLOS)
[Ntx, Nrx, Npath] = size(rayInfo.paths);

BistaticConstant = 8.734510731023729e-08;
LOSConstant = 2.744027474528536e-06;
Beta = 51.3127;

h = zeros(Ntx, Nrx, Npath, 2);
h0 = zeros(Npath, 2);

for i = 1:Ntx
    for j = 1:Nrx
        for p = 1:Npath
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
        h(i, j, :, :) = sortrows(h0, 2);
    end
end
h = squeeze(h);
% [Ntx, Nrx, Npath, 2] = size(h);

R = real(h(:, 2));
A = h(:, 1);
