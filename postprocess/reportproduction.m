function [VPI,oilrecovery,cumulateoil,watercut,countime]=reportproduction...
         (countime,VPI,wells,f_elem,cont,oilrecovery,cumulateoil,watercut,...
         q,dt,porousarea,S_old, bedge, bflux)
    global bcflag

if any(bcflag(:,1) > 201) % 201 flags for no-flow bc
    % modified by JCT
    % production faces
    bface  = bflux > 0;  
    fprod  = bflux(bface);
    pcells = bedge(bface,3);
    % injectin faces
    bface  = bflux < 0;  
    finj  = abs(bflux(bface));
    % pore volume
    pv = sum(porousarea);
    % compute variables
    sumpvi    = dt * sum( finj ) / pv ;
    sumcumoil = dt * sum( abs(f_elem(pcells) - 1) .* fprod );
    sumvol    = sum(porousarea(pcells));
    sumwater  = porousarea(pcells)*S_old(pcells);
else   
    sumpvi = 0;
    sumcumoil = 0;
    sumvol  = 0;
    sumwater=0;

    for iwell = 1:size(wells,1)
        
        if wells(iwell,3) == 0 % no poço produtor
            
            sumvol = sumvol + porousarea(wells(iwell,1));

            sumwater = sumwater + S_old(wells(iwell,1))*porousarea(wells(iwell,1));
                        
            f_o = 1 - f_elem(wells(iwell,1));
            
            sumcumoil = sumcumoil - f_o*q(wells(iwell,1))*dt;

        elseif wells(iwell,3) ~= 0 % no poço injetor
            
            sumpvi = sumpvi + f_elem(wells(iwell,1))*q(wells(iwell,1))*dt*(1/sum(porousarea));

        end
        
    end
end
% modified by JCT
countime(cont+1)    = countime(cont) + dt;
VPI(cont+1)         = VPI(cont) + sumpvi;
cumulateoil(cont+1) = cumulateoil(cont) + sumcumoil;
oilrecovery(cont+1) = 1 - (sumwater/sumvol);
watercut(cont+1)    = sumwater/sumvol;

end
