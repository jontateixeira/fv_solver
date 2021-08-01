function [influx,bflux]=flowratecalc(p,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,mobility)

global coord esurn1 esurn2 bedge inedge centelem bcflag af physmodel elem

% pre-aloca��o dos vetores
influx=zeros(size(inedge,1),1); % vetor vaz�o na face interior
bflux=zeros(size(bedge,1),1); % vetor vaz�o na face de contorno

nd=sum(af==0);
if strcmp(physmodel,'d')==1
    cg=2; 
else
    cg=1;
end

for iface=1:size(inedge,1)
           
    plef=p(inedge(iface,3)); %indice do elemento a direita da aresta i
    prel=p(inedge(iface,4)); %indice do elemento a esquerda da aresta i
    % interpolando os n�s (ou v�rtices) 
    nec1=esurn2(inedge(iface,1)+1)-esurn2(inedge(iface,1));
    p1=0;
    % calculo da press�o no n� "inedge(iface,1)"
    if nflag(inedge(iface,1),1)>200
        if (nflag(inedge(iface,1),1)>201)&&(nflag(inedge(iface,1),1)<300)
            for j=1:nec1
                element1=esurn1(esurn2(inedge(iface,1))+j);
                p1=p1+w(esurn2(inedge(iface,1))+j)*p(element1);
            end
            p1=p1+s(inedge(iface,1),1);
        else
            for j=1:nec1
                element1=esurn1(esurn2(inedge(iface,1))+j);
                p1=p1+w(esurn2(inedge(iface,1))+j)*p(element1);
            end
        end
        
    else
        p1=nflag(inedge(iface,1),2);
    end
    
    % calculo da press�o no "inedge(i,2)"
    nec2=esurn2(inedge(iface,2)+1)- esurn2(inedge(iface,2));
    p2=0;
    if nflag(inedge(iface,2),1)>200
        if (nflag(inedge(iface,2),1)>201)&&(nflag(inedge(iface,2),1)<300)
            for j=1:nec2
                element2=esurn1(esurn2(inedge(iface,2))+j);
                p2=p2+w(esurn2(inedge(iface,2))+j)*p(element2);
            end
            p2=p2+s(inedge(iface,2),1);
        else
            for j=1:nec2
                element2=esurn1(esurn2(inedge(iface,2))+j);
                p2=p2+w(esurn2(inedge(iface,2))+j)*p(element2);
            end
            
        end
        
    else
        p2=nflag(inedge(iface,2),2);
    end 
    
    % O elemento a direita ser� o que estiver na posi��o indicada pelo
    % inedge(i,4), ou seja p(rel)=p(pmps(inedge(i,4)). Se a linha inedge(i,4)
    % tiver mais de uma press�o, significa que � uma fratura, ent�o eu
    % tenho que ser se o elemento a esquerda � uma fratur tamb�m ou n�o. Se
    % for, ent�o eu calculo o fluxo com as m�dias das press�es das duas. Se
    % n�o for, ent�o eu vejo qual a press�o da fratura que � mais pr�xima
    % da press�o do elemento � esquerda e calculo o fluxo com ela.
     
    % calculo das vaz�es
    influx(iface)=mobility(iface)*Kde(iface)*(prel-plef-Ded(iface)*(p2-p1));
        
end

for ifacont=1:size(bedge,1)
    
    lef=bedge(ifacont,3);
    C=centelem(lef,:); % baricentro do elemento a esuqerda
    nor=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    if nor==0
        nor=af(elem(bedge(ifacont,3),5))/cg;
    end
    if bedge(ifacont,5)<200 % se os n�s esteverem na fronteira de DIRICHLET
        if bedge(ifacont,1)~=bedge(ifacont,2)
            c1=nflag(bedge(ifacont,1),2);
            c2=nflag(bedge(ifacont,2),2);
            bflux( ifacont )=-(Kn(ifacont)/(Hesq(ifacont)*nor))*(((C-coord(bedge(ifacont,2),:)))*(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:))'*c1+...
                (C-coord(bedge(ifacont,1),:))*(coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:))'*c2-(nor^2)*p(lef))+(c2-c1)*Kt(ifacont);
            bflux( ifacont )=mobility( ifacont+size(inedge,1) )*bflux( ifacont );
        else
            c1=nflag(bedge(ifacont,1),2);
            bflux( ifacont )=-(Kn(ifacont)/(Hesq(ifacont)))*(c1-p(lef));
            bflux( ifacont )=mobility( ifacont+size(inedge,1) )*bflux( ifacont );
        end
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        bflux(ifacont)= nor*bcflag(r,2);
    end
    
end

end