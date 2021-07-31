function [y] = harmonicopoint(kmap)

global inedge bedge coord elem centelem region

y=zeros(size(inedge,1)+size(bedge,1),3);
Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
for ifacont=1:size(bedge,1)
    y(ifacont,:)= 0.5*( coord(bedge(ifacont,1),:)+ coord(bedge(ifacont,2),:)); 
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % calculo das projeções ortogonais sobre a face
    %[hrel]=projortogonal(rel,inedge(iface,2), inedge(iface,1));
    
    %[hlef]=projortogonal(lef,inedge(iface,1), inedge(iface,2));
    %Determinação das alturas dos centróides dos elementos
    
    %Do ponto do início da aresta até o centro da célula da direita
    vd2=centelem(rel,:)-coord(inedge(iface,1),:);
    cd=cross(vd1,vd2);
    hrel=norm(cd)/norm(vd1); % altura a direita
    
    %Do ponto do início da aresta até o centro da célula da direita
    ve2=centelem(lef,:)-coord(inedge(iface,1),:);
    ce=cross(vd1,ve2);
    hlef=norm(ce)/norm(vd1); % altura a esquerda
    
    % tensor do elemento a esquerda
    
    Klef(1,1)=kmap(region(lef),2);
    Klef(1,2)=kmap(region(lef),3);
    Klef(2,1)=kmap(region(lef),4);
    Klef(2,2)=kmap(region(lef),5);
    
    % tensor do elemento a direita
    
    Krel(1,1)=kmap(region(rel),2);
    Krel(1,2)=kmap(region(rel),3);
    Krel(2,1)=kmap(region(rel),4);
    Krel(2,2)=kmap(region(rel),5);
    
    % calculo das constantes normais em cada face interna
    Knlef=dot(R*vd1',Klef*(R*vd1')/norm(vd1)^2);
    
    Knrel=dot((R*(-vd1')),Krel*(R*(-vd1'))/norm(vd1)^2);
    % calculo dos pontos armonicos
    y(iface+size(bedge,1),:)=(hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)'+...
        hlef*hrel*(Klef'-Krel')*(R*(vd1/norm(vd1))'))/(hrel*Knlef+hlef*Knrel);
    
%=========================================================================%
%     calculo dos pontos que saim fora da face interna nos pontos
%     intermeadiarios da discontinuidade, sobre na discontinuidade
%     horizontal.
%     if strcmp(benchmark,'benchmar5_6') || strcmp(benchmark,'benchmar5_7')
%         r1= find(N(inedge(iface,1),:)~=0 & N(inedge(iface,1),:)~=iface);
%         a1=N(inedge(iface,1),r1);
%         
%         s1= find(N(inedge(iface,2),:)~=0 & N(inedge(iface,2),:)~=iface );
%         b1=N(inedge(iface,2),s1);
%         
%         vetortotal= [a1 b1];
%         
%         for j=vetortotal
%             
%             if (j<size(inedge,1) || j==size(inedge,1))
%                 a1=inedge(j,1);
%                 a2=inedge(j,2);
%                 if (((coord(a1,1)< y(iface+size(bedge,1),1) || abs(coord(a1,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a2,1) || abs(y(iface+size(bedge,1),1)- coord(a2,1))<1e-10)) &&...
%                         ((coord(a1,2)< y(iface+size(bedge,1),2)|| abs(coord(a1,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a2,2) || abs(y(iface+size(bedge,1),2)-coord(a2,2))<1e-10))) ||...
%                         (((coord(a2,1)< y(iface+size(bedge,1),1) || abs(coord(a2,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a1,1) || abs(y(iface+size(bedge,1),1)- coord(a1,1))<1e-10)) &&...
%                         ((coord(a2,2)< y(iface+size(bedge,1),2)|| abs(coord(a2,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a1,2) || abs(y(iface+size(bedge,1),2)-coord(a1,2))<1e-10)))
%                     %                levando ao ponto interior
%                     y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
%                     break
%                 end
%             else
%                 a1=bedge(j-size(inedge,1),1);
%                 a2=bedge(j-size(inedge,1),2);
%                 if (((coord(a1,1)< y(iface+size(bedge,1),1) || abs(coord(a1,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a2,1) || abs(y(iface+size(bedge,1),1)- coord(a2,1))<1e-10)) &&...
%                         ((coord(a1,2)< y(iface+size(bedge,1),2)|| abs(coord(a1,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a2,2) || abs(y(iface+size(bedge,1),2)-coord(a2,2))<1e-10))) ||...
%                         (((coord(a2,1)< y(iface+size(bedge,1),1) || abs(coord(a2,1)- y(iface+size(bedge,1),1))<1e-10)  && (y(iface+size(bedge,1),1)< coord(a1,1) || abs(y(iface+size(bedge,1),1)- coord(a1,1))<1e-10)) &&...
%                         ((coord(a2,2)< y(iface+size(bedge,1),2)|| abs(coord(a2,2)- y(iface+size(bedge,1),2))<1e-10) && (y(iface+size(bedge,1),2) < coord(a1,2) || abs(y(iface+size(bedge,1),2)-coord(a1,2))<1e-10)))
%                     
%                     %               levando ao ponto interior
%                     y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
%                     
%                     break
%                 end
%             end
%         end
%         
%     elseif strcmp(benchmark,'edqueiroz')
%         a=0.5*(coord(inedge(iface,1),:)+coord(inedge(iface,2),:)); % ponto medio da face
%         
%         if norm(a-coord(inedge(iface,1),:))<norm(a-y(iface+size(bedge,1),:))
%             y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
%             %y(iface+size(bedge,1),:)= a;
%         end
%      else
%         % ative este calculo em malhas severamente distorcidas ja que
%          a=0.5*(coord(inedge(iface,1),:)+coord(inedge(iface,2),:)); % ponto medio da face
%          % este calculo eh usad por varios autores na literatura
%          %a= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel);
%          
%          if norm(a-coord(inedge(iface,1),:))<norm(a-y(iface+size(bedge,1),:))
%            %y(iface+size(bedge,1),:)= (hrel*Knlef*centelem(lef,:)'+ hlef*Knrel*centelem(rel,:)')/(hrel*Knlef+hlef*Knrel); % ativar para obter a convergencia Go e Wu 2010
%            y(iface+size(bedge,1),:)= a;
%          end
%         
%            end
   %======================================================================% 
end

end