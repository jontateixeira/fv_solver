function [S_old] = firstorderstandard(S_old,influx,bflux,q,f_elem,dt,wells,...
                                      S_cont,nw,no,elem,inedge,bedge,pormap,...
                                      elemarea,satlimit,visc,region)
global problemData
iter = 1;
if problemData.timestep.dvtol < -1
    iter = problemData.timestep.mtime(1);
    fprintf('* Perform %u substeps in saturation equation\n', iter);
end
for s=1:iter
RHS = zeros(size(elem,1),1);

for iface = 1 : size(inedge,1) 
    
    lef = inedge(iface,3); % elemento a esquerda
    rel = inedge(iface,4); % elemento a direita
    
    ve_mais = (influx(iface) + abs(influx(iface)))/2;
    ve_menos = (influx(iface) - abs(influx(iface)))/2;
    
    RHS(rel) = RHS(rel) + ve_mais*f_elem(lef) + ve_menos*f_elem(rel);
    RHS(lef)  = RHS(lef)  - ve_mais*f_elem(lef) - ve_menos*f_elem(rel);
    
end
% calculo dos contribuições nos poços 
if max(wells)~=0
    for iw=1:size(wells,1)
        prodwell=wells(iw,3);
        if prodwell==0
            sink = f_elem(wells(iw,1))*q(wells(iw,1));
            RHS(wells(iw,1)) = RHS(wells(iw,1)) + sink;
        end
    end
% calculo dos contribuições nos contornos onde tem fluxo prescrito
% diferente de zero
else
    for ifacont = 1:size(bedge,1)
        lef=bedge(ifacont,3);
        if (bedge(ifacont,5)>201)&&(bedge(ifacont,5)<300)
            % calculo do fluxo fracional no contorno com saturação
            % prescrito diferente de zero.
            if bflux(ifacont)<0
                Krw2 = ((S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^nw;
                Kro2 = ((1 - S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^no;
                L2 = Krw2/visc(1) + Kro2/visc(2);
                f_cont = (Krw2/visc(1))/L2;
                RHS(lef) = RHS(lef) - f_cont * bflux(ifacont);
            else
                RHS(lef) = RHS(lef) - f_elem(bedge(ifacont,3)) * bflux(ifacont);
            end
            
        elseif bedge(ifacont,5)<200
            
            if bflux(ifacont)<0
                Krw2 = ((S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^nw;
                Kro2 = ((1 - S_cont - satlimit(1))/(1-satlimit(1)-satlimit(2)))^no;
                L2 = Krw2/visc(1) + Kro2/visc(2);
                f_cont = (Krw2/visc(1))/L2;
                RHS(lef) = RHS(lef) - f_cont * bflux(ifacont);
            else
                RHS(lef) = RHS(lef) - f_elem(bedge(ifacont,3)) * bflux(ifacont);
            end
            
        end
    end
end

for i = 1:size(S_old,1)
    por=pormap(region(i));
    if S_old(i)~=S_cont
        S_old(i,1) = S_old(i,1) + (dt*RHS(i))/(por*elemarea(i));
    end
end

end