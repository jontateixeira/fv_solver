function [coefficient]=coefficientLPSangle(kmap,nflag)
global inedge bedge coord elem centelem eps region

Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0; 0 0 0];

for ifacont=1:size(bedge,1)
    
    lef=bedge(ifacont,3);
    klef=elementype(lef);
    
    IJ=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    mid=0.5*(coord(bedge(ifacont,2),:)+coord(bedge(ifacont,1),:));
    clef=centelem(lef,:);
    normIJ=norm(IJ);
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(region(bedge(ifacont,3)),2);
    Klef(1,2)=kmap(region(bedge(ifacont,3)),3);
    Klef(2,1)=kmap(region(bedge(ifacont,3)),4);
    Klef(2,2)=kmap(region(bedge(ifacont,3)),5);
    
    %% ej(iface) do elemento a esquerda
    
    ve2= Klef'*(R*(IJ/normIJ)');
    
%     ve2=mid-clef; ve2=ve2';
    
    ve21=ve2/norm(ve2);
    
    %% percorrendo todos as faces dos elemento "lef"
    auxvetor=elem(lef,1:klef);
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            no1=i;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(lef,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve21)/(norm(vj)*norm(ve21))>1 && abs(1-dot(vj,ve21)/(norm(vj)*norm(ve21)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve21)/(norm(vj)*norm(ve21)));
            end
            vi=coord(no1,:)-centelem(lef,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve21);
            auxquadrant2= cross(ve21,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve21)/(norm(vi)*norm(ve21))>1 && abs(1-dot(vi,ve21)/(norm(vi)*norm(ve21)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve21)/(norm(vi)*norm(ve21)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,ifacont)=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,ifacont)=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,ifacont)=aux11;
                coefficient(1,4,ifacont)=aux12;
                
                % verificando a identidade
                coefficient(1,5,ifacont)=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
                if abs(coefficient(1,5,ifacont))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        else
            no1=i;
            no2=auxvetor(j+1);
            vj= coord(no2,:)-centelem(lef,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve21)/(norm(vj)*norm(ve21))>1 && abs(1-dot(vj,ve21)/(norm(vj)*norm(ve21)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve21)/(norm(vj)*norm(ve21)));
            end
            vi= coord(no1,:)-centelem(lef,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,ve21);
            auxquadrant2= cross(ve21,vj);
            % evita que aparição de numeros complexos
            if dot(vi,ve21)/(norm(vi)*norm(ve21))>1 && abs(1-dot(vi,ve21)/(norm(vi)*norm(ve21)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve21)/(norm(vi)*norm(ve21)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,ifacont)=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,ifacont)=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,ifacont)=aux11;
                coefficient(1,4,ifacont)=aux12;
                
                % verificando a identidade
                coefficient(1,5,ifacont)=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+(coord(aux12,:)-...
                    centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
                if abs(coefficient(1,5,ifacont))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
            end
            
        end
        j=j+1;
    end
    
    clear aux12 aux11
end

%% Faces interiores
for iface=1:size(inedge,1)
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    klef=elementype(lef);
    krel=elementype(rel);
    
    IJ=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    normIJ=norm(coord(inedge(iface,2),:)-coord(inedge(iface,1),:));
    LR=centelem(lef,:)-centelem(rel,:);
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(region(inedge(iface,3)),2);
    Klef(1,2)=kmap(region(inedge(iface,3)),3);
    Klef(2,1)=kmap(region(inedge(iface,3)),4);
    Klef(2,2)=kmap(region(inedge(iface,3)),5);
    
    %tensor a elemento a direita
    
    Krel(1,1)=kmap(region(inedge(iface,4)),2);
    Krel(1,2)=kmap(region(inedge(iface,4)),3);
    Krel(2,1)=kmap(region(inedge(iface,4)),4);
    Krel(2,2)=kmap(region(inedge(iface,4)),5);
    
    %% ej(iface) do elemento a esquerda
    
    % terceira maneira de calcular o kn
%     if nflag(inedge(iface,1),1)<300 | nflag(inedge(iface,2),1)<300
%         ve2=LR; ve2=ve2';
%     else
        ve2=Klef*(R*(IJ/normIJ)');
%     end
    
    ve21=ve2/norm(ve2);
    % faces na que compoem o elemento
    auxvetor=elem(lef,1:klef);
    
    j=1;
    for i=auxvetor
        if j==length(auxvetor)
            no1=i;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(lef,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve21)/(norm(vj)*norm(ve21))>1 && abs(1-dot(vj,ve21)/(norm(vj)*norm(ve21)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve21)/(norm(vj)*norm(ve21)));
            end
            vi=coord(no1,:)-centelem(lef,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve21);
            auxquadrant2= cross(ve21,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve21)/(norm(vi)*norm(ve21))>1 && abs(1-dot(vi,ve21)/(norm(vi)*norm(ve21)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve21)/(norm(vi)*norm(ve21)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,iface+size(bedge,1))=aux11;
                coefficient(1,4,iface+size(bedge,1))=aux12;
                
                % verificando a identidade
                coefficient(1,5,iface+size(bedge,1))=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))- ve2');
                if abs(coefficient(1,5,iface+size(bedge,1)))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        else
            no1=i;
            no2=auxvetor(j+1);
            vj= coord(no2,:)-centelem(lef,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,ve21)/(norm(vj)*norm(ve21))>1 && abs(1-dot(vj,ve21)/(norm(vj)*norm(ve21)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vj,ve21)/(norm(vj)*norm(ve21)));
            end
            vi= coord(no1,:)-centelem(lef,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vj,ve21);
            auxquadrant2= cross(ve21,vi);
            % evita que aparição de numeros complexos
            if dot(vi,ve21)/(norm(vi)*norm(ve21))>1 && abs(1-dot(vi,ve21)/(norm(vi)*norm(ve21)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vi,ve21)/(norm(vi)*norm(ve21)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                
                aux11=no1;
                aux12=no2;
                % atribuindo valores a os coeficientes
                coefficient(1,1,iface+size(bedge,1))=(norm(ve2))*sin(thetalef2)/(norm(coord(aux11,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                coefficient(1,2,iface+size(bedge,1))=(norm(ve2))*sin(thetalef1)/(norm(coord(aux12,:)-centelem(lef,:))*sin(thetalef1+thetalef2));
                
                % indexando as faces respetivamente
                coefficient(1,3,iface+size(bedge,1))=aux11;
                coefficient(1,4,iface+size(bedge,1))=aux12;
                
                % verificando a identidade
                coefficient(1,5,iface+size(bedge,1))=norm((coord(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
                    (coord(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))- ve2');
                if abs(coefficient(1,5,iface+size(bedge,1)))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
            end
            
        end
        j=j+1;
    end
    
    clear aux12 aux11
    
    %% Elemento a direita
%     if nflag(inedge(iface,1),1)<300 | nflag(inedge(iface,2),1)<300
%         vd2=-LR; vd2=vd2';
%     else
        vd2=Krel*(R*(-IJ/normIJ)');
%     end
    auxvetor=elem(rel,1:krel);
    vd21=vd2/norm(vd2);
    j=1;
    for ii=auxvetor
        if j==length(auxvetor)
            no1=ii;
            no2=auxvetor(1);
            vj=coord(no2,:)-centelem(rel,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,vd21)/(norm(vj)*norm(vd21))>1 && abs(1-dot(vj,vd21)/(norm(vj)*norm(vd21)))<eps
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vj,vd21)/(norm(vj)*norm(vd21)));
            end
            vi=coord(no1,:)-centelem(rel,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,vd21);
            auxquadrant2= cross(vd21,vj);
            % evita que aparição de numeros complexos
            if dot(vi,vd21)/(norm(vi)*norm(vd21))>1 && abs(1-dot(vi,vd21)/(norm(vi)*norm(vd21)))<eps
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vi,vd21)/(norm(vi)*norm(vd21)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                
                auxi=no1;
                auxj=no2;
                % atribuindo valores a os coeficientes
                coefficient(2,1,iface+size(bedge,1))=(norm(vd2))*sin(thetarel2)/(norm(coord(auxi,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                coefficient(2,2,iface+size(bedge,1))=(norm(vd2))*sin(thetarel1)/(norm(coord(auxj,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                
                % indexando as faces respetivamente
                coefficient(2,3,iface+size(bedge,1))=auxi;
                coefficient(2,4,iface+size(bedge,1))=auxj;
                
                % verificando a identidade
                coefficient(2,5,iface+size(bedge,1))=norm((coord(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
                    (coord(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))- vd2');
                if abs(coefficient(2,5,iface+size(bedge,1)))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                    
                end
                
            end
            
        else
            no1=ii;
            no2=auxvetor(j+1);
            vj=coord(no2,:)-centelem(rel,:);
            vj=vj/norm(vj);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vj,vd21)/(norm(vj)*norm(vd21))>1 && abs(1-dot(vj,vd21)/(norm(vj)*norm(vd21)))<eps
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vj,vd21)/(norm(vj)*norm(vd21)));
            end
            
            vi= coord(no1,:)-centelem(rel,:);
            vi=vi/norm(vi);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vi,vd21);
            auxquadrant2= cross(vd21,vj);
            % evita que aparição de numeros complexos
            if dot(vi,vd21)/(norm(vi)*norm(vd21))>1 && abs(1-dot(vi,vd21)/(norm(vi)*norm(vd21)))<eps
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vi,vd21)/(norm(vi)*norm(vd21)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps ||...
                    abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                
                auxi=no1;
                auxj=no2;
                % atribuindo valores a os coeficientes
                coefficient(2,1,iface+size(bedge,1))=(norm(vd2))*sin(thetarel2)/(norm(coord(auxi,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                coefficient(2,2,iface+size(bedge,1))=(norm(vd2))*sin(thetarel1)/(norm(coord(auxj,:)-centelem(rel,:))*sin(thetarel1+thetarel2));
                
                % indexando as faces respetivamente
                coefficient(2,3,iface+size(bedge,1))=auxi;
                coefficient(2,4,iface+size(bedge,1))=auxj;
                
                % verificando a identidade
                coefficient(2,5,iface+size(bedge,1))=norm((coord(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
                    (coord(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))- vd2');
                if abs(coefficient(2,5,iface+size(bedge,1)))<eps
                else
                    
                    disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
                    
                end
                
            end
        end
        j=j+1;
        
    end
    
    clear auxj auxi
    
end

end