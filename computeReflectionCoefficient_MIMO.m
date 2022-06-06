function Gm = computeReflectionCoefficient_MIMO(app, wall, d)
            % Reflection coefficient depending on wall properties and
            % frequency
            res = d;
            G = @(Z1, Z2, cosi, cost) (Z2*cosi-Z1*cost)./(Z2*cosi+Z1*cost); % transmission coefficient for an infinite wall
            sini = @(cosi) sqrt(1-cosi^2);
            sint = @(e1, e2, sini) (sini*sqrt(e1/e2));
            cost = @(sint) sqrt(1-sint^2);
            width = @(l, cost) l/cost;
            Gm = @(G, gamma_m, beta, s, sint, sini, f) (G - (1-G^2) .* G.*exp(-2*gamma_m*s).*exp(1i*beta*2*s.*sint.*sini) ...
                ./(1-G.^2.*exp(-2*gamma_m*s).*exp(1i*beta*2*s.*sint.*sini))); % reflection coefficient for a wall

            % Steering vector
            %ci = abs(res * [wall.ny; -wall.nx] / norm(res));
            ci = abs([d(1,1,1,1), d(1,1,1,2)] * [wall.ny; -wall.nx] / norm([d(1,1,1,1), d(1,1,1,2)]));
            si = sini(ci);
            st = sint(1, wall.Permittivity, si);
            ct = cost(st);
            Gm = Gm(G(app.Z0, wall.Z, ci, ct), wall.gamma, ...
                app.beta, width(wall.Width, ct), st, si, app.f);
        end