function [S_old] = firstorderstdseqimp(S_old,influx,bflux,dt,wells,...
                          q,nw,no,elem,bedge,inedge,elemarea,pormap,region)
         
C = sparse(size(elem,1),size(elem,1));

t1=1;t2=1;
if isempty(wells)==0
    for i=1:size(wells,1)
       if wells(i,3)>300 && wells(i,3)<400
           injecelem(1,t1)=wells(i,1);
           t1=t1+1;
       elseif wells(i,3)==0
           producelem(1,t2)=wells(i,1);
           t2=t2+1;
       end
    end
else
    injecelem = bedge(bflux<0,3); S_old(injecelem) = 1 - eps;
    producelem = bedge(bflux>0,3);
end

%--------------------------------------------------------------------------
%Boundary edges
for i = 1:size(bedge,1)
    C(bedge(i,3),bedge(i,3)) =  C(bedge(i,3),bedge(i,3)) + ...
        bflux(i)/(elemarea(bedge(i,3))*pormap(region(bedge(i,3))));
end

%--------------------------------------------------------------------------
%Internal edges
for i = 1:size(inedge,1)    
    if influx(i) >= 0
        C(inedge(i,3),inedge(i,3)) =  C(inedge(i,3),inedge(i,3)) - ...
            influx(i)/(elemarea(inedge(i,3))*pormap(region(inedge(i,3))));            
        C(inedge(i,4),inedge(i,3)) =  C(inedge(i,4),inedge(i,3)) + ...
            influx(i)/(elemarea(inedge(i,4))*pormap(region(inedge(i,4))));
    else
        C(inedge(i,3),inedge(i,4)) =  C(inedge(i,3),inedge(i,4)) - ...
            influx(i)/(elemarea(inedge(i,3))*pormap(region(inedge(i,3))));
        C(inedge(i,4),inedge(i,4)) =  C(inedge(i,4),inedge(i,4)) + ...
            influx(i)/(elemarea(inedge(i,4))*pormap(region(inedge(i,4))));   
    end    
end 

%Calculo termo fonte explicito
for i=1:size(producelem,2)
    Qprod = q(producelem(i))/(elemarea(producelem(i))*pormap(region(producelem(i))));
    C(producelem(i),producelem(i)) = C(producelem(i),producelem(i)) + Qprod;
end
for i=1:size(injecelem,2)
    Qinj = q(injecelem(i))/(elemarea(injecelem(i))*pormap(region(injecelem(i))));
    C(injecelem(i),injecelem(i)) = C(injecelem(i),injecelem(i)) + Qinj;
end

S0 = S_old;
dt0 = dt;
d_t = 0;
pdt=0;

while d_t < dt0
    
    fprintf('%u%%',pdt);
    ref=0; 
    while ref == 0
                
        [ Sw_n1k1, r ] = NewtonRaphson( S0, nw, no, dt, C );

        if (max(Sw_n1k1) > 1)||(min(Sw_n1k1) < 0)||(max(abs(r)) > 10^-3)

            dt = dt/2;
            ref = 0;
                       
        else

            ref = ref + 1;
            d_t = d_t + dt;
%             dt = dt0 - d_t;

            if dt > dt0 - d_t
                dt = dt0 - d_t;
            end

        end

    end
    
    S0 = Sw_n1k1;
    
    if pdt<10
        fprintf('\b\b');
    else
        fprintf('\b\b\b');
    end
    pdt = floor(100*(d_t/dt0));
        
end

S_old = S0;
fprintf('%u%%\n',pdt);

end

