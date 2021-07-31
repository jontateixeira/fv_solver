function [Sw_n1k1] = firstorderstdfrat(Sw_n1k1,S0,influx,f_elem,dt,dt_frat,elem,inedge)

global af porousarea

porsel = porousarea';

RHS=zeros(size(elem,1),1);

nd=sum(af==0);

t=1;
for i=1:size(inedge,1)
   if elem(inedge(i,3),5)>nd || elem(inedge(i,4),5)>nd
      ve(t,1)=i;
      t=t+1;
   end
end

z=find(elem(:,5)>nd);

Sw_n1k1(z)=S0(z);

an=(round(dt/dt_frat)+1);
if an>10
    n=an;
else
    n=10;
end
dt1=dt/n;

for j=1:n

    for i=1:size(ve,1)
        
        f=ve(i);

        lef = inedge(f,3); % elemento a esquerda
        rel = inedge(f,4); % elemento a direita

        ve_mais = (influx(f) + abs(influx(f)))/2;
        ve_menos = (influx(f) - abs(influx(f)))/2;

        RHS(rel) = RHS(rel) + ve_mais*f_elem(rel) + ve_menos*f_elem(rel);
        RHS(lef)  = RHS(lef)  - ve_mais*f_elem(lef) - ve_menos*f_elem(lef);

    end
    
    Sw_n1k1 = Sw_n1k1 + dt1*(RHS./porsel);

end

end