function [M,I]=assemblematrix_NLFV(pinterp,csi,wells,mobility)

global inedge coord bedge bcflag elem elemarea source

I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;

%%
sumvol=0;
for iw=1:size(wells,1)      
    if wells(iw,5) == 0 % Caso haja fluxo prescrito em algum poço.
        I(wells(iw,1)) = wells(iw,6)*elemarea(wells(iw,1));
        sumvol = sumvol + elemarea(wells(iw,1));
    end
end
if sumvol>0
    I=I./sumvol;
end

I=I+source;

for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);

    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));

    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)-normcont*bcflag(r,2);
    else
        %% calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)

        silef = mobility(ifacont)*normcont*(csi(1,1,ifacont)*pinterp(csi(1,3,ifacont))+...
            csi(1,2,ifacont)*pinterp(csi(1,4,ifacont)));

        Alef = mobility(ifacont)*normcont*(csi(1,1,ifacont)+csi(1,2,ifacont));

        %% implementação da matriz global no contorno
        M(lef,lef) = M(lef,lef) + Alef;
        I(lef,1) = I(lef,1) + silef;
    end
end

% I(wells(iw,1))= 1*elemarea(wells(iw,1)); !!!!!!!!!!!!!!!!!!!!!!!!!!!      
% injeta um m3 de agua por dia (d) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma= sqrt(vd1(1,1)^2+vd1(1,2)^2);
    ifactual=iface+size(bedge,1);
    
    % calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    silef=csi(1,1,ifactual)*pinterp(csi(1,3,ifactual))+...
        csi(1,2,ifactual)*pinterp(csi(1,4,ifactual));
    
    % direita
    sirel= csi(2,1,ifactual)*pinterp(csi(2,3,ifactual))+...
        csi(2,2,ifactual)*pinterp(csi(2,4,ifactual));
      
    wlef=(abs(sirel)+1e-16)/(abs(silef)+abs(sirel)+2*1e-16);
    wrel=(abs(silef)+1e-16)/(abs(silef)+abs(sirel)+2*1e-16);
    
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*wlef*(csi(1,1,ifactual)+csi(1,2,ifactual));
    ARR=norma*wrel*(csi(2,1,ifactual)+csi(2,2,ifactual));
    
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ mobility(ifactual)*ALL;
    M(lef,rel)=M(lef,rel)- mobility(ifactual)*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ mobility(ifactual)*ARR;
    M(rel,lef)=M(rel,lef)- mobility(ifactual)*ALL;    
end

for iw=1:size(wells,1)      
    if (wells(iw,5)>400)&&(wells(iw,5)<600) % Caso seja prescrita alguma pressão.
        M(wells(iw,1),:)=0*M(wells(iw,1),:);
        M(wells(iw,1),wells(iw,1))=1;
        I(wells(iw,1))=wells(iw,6);
    end
end

end