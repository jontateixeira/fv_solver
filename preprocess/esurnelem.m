function [esuel, esuel_coord,A,bound] = esurnelem(inedge)
global elem bedge coord centelem
%

proj=@(u,v)(dot(u,v)/dot(v,v))*v;
% a ordenação esta tudo em sentido antihorario 
esuel=zeros(size(elem,1),6);
bound=zeros(size(elem,1),1);
for i=1:size(elem,1)
   t=1;
   for j=1:size(inedge,1)
       if inedge(j,3)==i
           esuel(i,t)=inedge(j,4);
           t=t+1;
       elseif inedge(j,4)==i
           esuel(i,t)=inedge(j,3);
           t=t+1;
       end
   end
   for j=size(bedge,1)
      if bedge(j,3)==i
          bound(i)=j;
      end
   end
end    
% elemento repitidos da matriz 'esuel' somente ocorre quando o elemento
% forma parte do contorno.


% calculando as coordenadas dos elementos fantasma 
for i = 1:size(esuel,1)
    
    s=1;
    for kk=1:size(esuel,2) % loop sobre as colunas da matriz esuel
        
        if esuel(i,kk)~=i & esuel(i,kk)~=0 % evita o indice de um vetor seja zero
            % e excluye as proprierades geometrica do elemento atual 'i'
            
            esuel_coord(kk,:,i)=centelem(esuel(i,kk),:); % Ordena os centroides 
            % dos elementos na vizinhança do elemento em questão
            
            A(kk,:,i)= centelem(esuel(i,kk),:)-centelem(i,:); % diferença de centroides
            
        elseif esuel(i,kk)==i & esuel(i,kk)~=0 % evita o indice de um vetor seja zero
            % e incluye as proprierades geometrica do elemento atual 'i',
            % isto, acontece somente para elemento do contorno.
            for iface =s:size(bound,2)
                if bound(i,iface)~=0
                    
                 % esta pequena rutina calcula o centroide do elemento
                 % fantasma
                coordv1 = coord(bedge(bound(i,iface),1),:);
                coordv2 = coord(bedge(bound(i,iface),2),:);
                xL = centelem(i,:);
                
                p = proj((xL-coordv1),(coordv2-coordv1));
                aux3 = xL-(coordv1+p);
                
                esuel_coord(kk,:,i) = (xL-2*aux3)'; % centroide do elemento
                % fantasma
               
                A(kk,:,i)=(xL-2*aux3)-centelem(i,:); % diferença do centroido
                % do elemento 'i' a centroide do elemento fantasma
                s=s+1;
                break
                
                end
            end
        end
        
    end
    
end

end