function [coefficient,auxface]=coefficientPPSharmonicpoint(F,y,kmap)

global inedge bedge coord elem centelem eps region

Klef=zeros(3,3);
Krel=zeros(3,3);
R=[0 1 0; -1 0 0;0 0 0];
auxface=zeros(size(bedge,1)+size(inedge,1),4);
% eps=1e-5;

for ifacont=1:size(bedge,1)
    
    lef=bedge(ifacont,3);
    klef=elementype(lef);
    
    IJ=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normIJ=norm(coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:));
    %Essa é UMA maneira de construir os tensores.
    
    Klef(1,1)=kmap(region(bedge(ifacont,3)),2);
    Klef(1,2)=kmap(region(bedge(ifacont,3)),3);
    Klef(2,1)=kmap(region(bedge(ifacont,3)),4);
    Klef(2,2)=kmap(region(bedge(ifacont,3)),5);
    
    %% ej(iface) do elemento a esquerda
    
    ve2= Klef'*(R*(IJ/normIJ)');
    %% percorrendo todos as faces dos elemento "lef"
    auxvetor=F(lef,1:klef);
    j=1;
    
    ksii=1e30;
    ksij=1e30;
    for i=auxvetor
        if i~=auxvetor(length(auxvetor))
            vej=y(auxvetor(j+1),:)-centelem(lef,:);
            
            if dot(vej,ve2)/(norm(vej)*norm(ve2))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            end
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vei,ve2);
            % Estes condições evitam que o acos seja numero complexo.
            auxquadrant2= cross(ve2,vej);
            % evita que aparição de numeros complexos
            if dot(vei,ve2)/(norm(vei)*norm(ve2))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))|| ...
                (sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) ...
                && ((thetalef2 + thetalef1)<pi)
                ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                ksij=dot(cross(vei,ve2),[0,0,1])/dot(cross(vei,vej),[0 0 1]);
                aux11=i;
                aux12=auxvetor(j+1);
                            
            end
        else
            vej=y(auxvetor(1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if dot(vej,ve2)/(norm(vej)*norm(ve2))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            end
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vei,ve2);
            auxquadrant2= cross(ve2,vej);
            % evita que aparição de numeros complexos
            if dot(vei,ve2)/(norm(vei)*norm(ve2))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<eps
                % calculo do theta1
                thetalef1=acos(1);
            else
                thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))|| ...
                (sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) ...
                && ((thetalef2 + thetalef1)<pi)
                ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                ksij=dot(cross(vei,ve2),[0,0,1])/dot(cross(vei,vej),[0 0 1]);
                aux11=i;
                aux12=auxvetor(1);
            
            end
            
        end
        j=j+1;
    end
    
    if (ksii==1e30 && ksij==1e30) || (ksii>1e10 && ksij>1e10)
        ppp=1
        [ksii,ksij,aux11,aux12,auxy]=aroundfacelement(F,y,lef,ve2,klef,kmap);
        % atribuindo valores a os coeficientes
        coefficient(1,1,ifacont)=ksii;
        coefficient(1,2,ifacont)=ksij;
        
        % indexando as faces respetivamente
        coefficient(1,3,ifacont)=aux11;
        coefficient(1,4,ifacont)=aux12;
        
        % verificando a identidade
        coefficient(1,5,ifacont)=norm((auxy(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+ (auxy(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
        if abs(coefficient(1,5,ifacont))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
        end
    else
        % atribuindo valores a os coeficientes
        coefficient(1,1,ifacont)=ksii;
        coefficient(1,2,ifacont)=ksij;
        
        % indexando as faces respetivamente
        coefficient(1,3,ifacont)=aux11;
        coefficient(1,4,ifacont)=aux12;
        
        % verificando a identidade
        coefficient(1,5,ifacont)=norm((y(aux11,:)-centelem(lef,:))*coefficient(1,1,ifacont)+ (y(aux12,:)-centelem(lef,:))*coefficient(1,2,ifacont)-ve2');
        if abs(coefficient(1,5,ifacont))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
        end
    end
    [auxface]=calfacelement(lef,ifacont,aux11,aux12,bedge,inedge,auxface,1,2);
    
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
    
    %% elemento a esquerda
    
    ve2=Klef*(R*(IJ/normIJ)');
    % faces na que compoem o elemento
    auxvetor=F(lef,1:klef);
    ksii=1e30;
    ksij=1e30;
    j=1;
    for i=auxvetor
        if i~=auxvetor(length(auxvetor))
            vej=y(auxvetor(j+1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if (dot(vej,ve2)/(norm(vej)*norm(ve2)))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            end
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vei,ve2);
            auxquadrant2= cross(ve2,vej);
            % evita que aparição de numeros complexos
            if (dot(vei,ve2)/(norm(vei)*norm(ve2)))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<eps
                % calculo do theta1
                thetalef1=acos(1);
                
            else
                thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            end
            %% era >eps, tenga cuidado!!!!!!!!
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                ksij=dot(cross(vei,ve2),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                aux11=i;
                aux12=auxvetor(j+1);
            
                
            end
        else
            vej=y(auxvetor(1),:)-centelem(lef,:);
            % Estes condições evitam que o acos seja numero complexo.
            if (dot(vej,ve2)/(norm(vej)*norm(ve2)))>1 && abs(1-dot(vej,ve2)/(norm(vej)*norm(ve2)))<eps
                thetalef2=acos(1);
            else
                thetalef2=acos( dot(vej,ve2)/(norm(vej)*norm(ve2)));
            end
            
            vei=y(i,:)-centelem(lef,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vei,ve2);
            auxquadrant2= cross(ve2,vej);
            % evita que aparição de numeros complexos
            if (dot(vei,ve2)/(norm(vei)*norm(ve2)))>1 && abs(1-dot(vei,ve2)/(norm(vei)*norm(ve2)))<eps
                
                % calculo do theta1
                thetalef1=acos(1);
                
            else
                thetalef1=acos(dot(vei,ve2)/(norm(vei)*norm(ve2)));
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0)) && ((thetalef2 + thetalef1)<pi)
                ksii=dot(cross(ve2,vej),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                ksij=dot(cross(vei,ve2),[0 0 1])/dot(cross(vei,vej),[0 0 1]);
                aux11=i;
                aux12=auxvetor(1);
           
            end
            
        end
        j=j+1;
    end
    if (ksii==1e30 && ksij==1e30) || (ksii>1e10 && ksij>1e10)
        
        [ksii,ksij,aux11,aux12,auxy]=aroundfacelement(F,y,lef,ve2,klef,kmap);
        % atribuindo valores a os coeficientes
        coefficient(1,1,iface+size(bedge,1))=ksii;
        coefficient(1,2,iface+size(bedge,1))=ksij;
        
        % indexando as faces respetivamente
        coefficient(1,3,iface+size(bedge,1))=aux11;
        coefficient(1,4,iface+size(bedge,1))=aux12;
        
        % verificando a identidade
        coefficient(1,5,iface+size(bedge,1))=norm((auxy(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
            (auxy(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))-ve2');
        if abs(coefficient(1,5,iface+size(bedge,1)))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
        end
    else
        
        % atribuindo valores a os coeficientes
        coefficient(1,1,iface+size(bedge,1))=ksii;
        coefficient(1,2,iface+size(bedge,1))=ksij;
        
        % indexando as faces respetivamente
        coefficient(1,3,iface+size(bedge,1))=aux11;
        coefficient(1,4,iface+size(bedge,1))=aux12;
        
        % verificando a identidade
        coefficient(1,5,iface+size(bedge,1))=norm((y(aux11,:)-centelem(lef,:))*coefficient(1,1,iface+size(bedge,1))+...
            (y(aux12,:)-centelem(lef,:))*coefficient(1,2,iface+size(bedge,1))-ve2');
        if abs(coefficient(1,5,iface+size(bedge,1)))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
        end
    end
    [auxface]=calfacelement(lef,iface+size(bedge,1),aux11,aux12,bedge,inedge,auxface,1,2);
    
    clear ksii ksij aux11 aux12
    
    %% Elemento a direita
    vetor12=Krel*(R*(-IJ/normIJ)');
    auxvetorel=F(rel,1:krel);
     ksii=1e30;
     ksij=1e30;
    j=1;
    for ii=auxvetorel
        if ii~=auxvetorel(length(auxvetorel))
            vetorj=y(auxvetorel(j+1),:)-centelem(rel,:);
            % Estes condições evitam que o acos seja numero complexo.
            if (dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))>1 && abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))<eps
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)));
            end
            
            vetori=y(ii,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetori,vetor12);
            auxquadrant2= cross(vetor12,vetorj);
            % evita que aparição de numeros complexos
            if (dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))>1 && abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))<eps
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                ksii=dot(cross(vetor12,vetorj),[0 0 1])/dot(cross(vetori,vetorj),[0 0 1]);
                ksij=dot(cross(vetori,vetor12),[0 0 1])/dot(cross(vetori,vetorj),[0 0 1]);
                auxi=ii;
                auxj=auxvetorel(j+1);
                          
            end
            
        else
            vetorj=y(auxvetorel(1),:)-centelem(rel,:);
            % Estes condições evitam que o acos seja numero complexo.
            if (dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))>1 && abs(1-dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)))<eps
                % calculo do theta2
                thetarel2=acos(1);
            else
                thetarel2=acos( dot(vetorj,vetor12)/(norm(vetorj)*norm(vetor12)));
            end
            vetori=y(ii,:)-centelem(rel,:);
            % analiza que o K.n pertece ao primeiro quadrante
            auxquadrant1= cross(vetori,vetor12);
            auxquadrant2= cross(vetor12,vetorj);
            % evita que aparição de numeros complexos
            if (dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))>1 && abs(1-dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)))<eps
                % calculo do theta1
                thetarel1=acos(1);
            else
                thetarel1=acos( dot(vetori,vetor12)/(norm(vetori)*norm(vetor12)));
                
            end
            if ((sign(auxquadrant1(1,3))==sign(auxquadrant2(1,3)) && (abs(auxquadrant1(1,3))>eps || abs(auxquadrant2(1,3))>eps))||(sign(auxquadrant1(1,3))==0 && sign(auxquadrant2(1,3))>0)||...
                    (sign(auxquadrant1(1,3))>0 && sign(auxquadrant2(1,3))==0))  && ((thetarel2 + thetarel1)<pi)
                ksii=dot(cross(vetor12,vetorj),[0 0 1])/dot(cross(vetori,vetorj),[0 0 1]);
                ksij=dot(cross(vetori,vetor12),[0 0 1])/dot(cross(vetori,vetorj),[0 0 1]);
                auxi=ii;
                auxj=auxvetorel(1);
           
            end
        end
        j=j+1;
        
    end
    if (ksii==1e30 && ksij==1e30)|| (ksii>1e10 && ksij>1e10)
        
        [ksii,ksij,auxi,auxj,auxy]=aroundfacelement(F,y,rel,vetor12,krel,kmap);
        % atribuindo valores a os coeficientes
        coefficient(2,1,iface+size(bedge,1))=ksii;
        coefficient(2,2,iface+size(bedge,1))=ksij;
        % indexando as faces respetivamente
        coefficient(2,3,iface+size(bedge,1))=auxi;
        coefficient(2,4,iface+size(bedge,1))=auxj;
        % verificando a identidade
        coefficient(2,5,iface+size(bedge,1))=norm((auxy(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
            (auxy(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))-vetor12');
        if abs(coefficient(2,5,iface+size(bedge,1)))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
            
        end
    else
        % atribuindo valores a os coeficientes
        coefficient(2,1,iface+size(bedge,1))=ksii;
        coefficient(2,2,iface+size(bedge,1))=ksij;
        % indexando as faces respetivamente
        coefficient(2,3,iface+size(bedge,1))=auxi;
        coefficient(2,4,iface+size(bedge,1))=auxj;
        % verificando a identidade
        coefficient(2,5,iface+size(bedge,1))=norm((y(auxi,:)-centelem(rel,:))*coefficient(2,1,iface+size(bedge,1))+...
            (y(auxj,:)-centelem(rel,:))*coefficient(2,2,iface+size(bedge,1))-vetor12');
        if abs(coefficient(2,5,iface+size(bedge,1)))<eps
        else
            
            disp('Não satisfaz a identidade (7) do artigo Gao and Wu (2013)');
            
        end
    end
    [auxface]=calfacelement(rel,iface+size(bedge,1),auxi,auxj,bedge,inedge,auxface,3,4);
    clear ksii ksij auxi auxj
    
end

end