function [M,I]=assemblematrix_MPFAQL(parameter,w,s,nflag,weightDMP,wells,mobility)

global inedge coord bedge bcflag elem elemarea esurn1 esurn2 source

I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;

%%
sumvol=0;
for iw=1:size(wells,1)      
    if wells(iw,1) ~= 0 && wells(iw,5) == 0 % Caso haja fluxo prescrito em algum poço.
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
        I(lef)=I(lef)- normcont*bcflag(r,2);
    else
        % implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ mobility(ifacont)*normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        nolef1=parameter(1,3,ifacont);
        nolef2=parameter(1,4,ifacont);
        % contribuições do nó 1
        if nflag(nolef1,1)<200
            I(lef,1)=I(lef,1)+ mobility(ifacont)*normcont*(parameter(1,1,ifacont)*nflag(nolef1,2));
        else

            for j=1:(esurn2(nolef1+1)-esurn2(nolef1))

                post_cont=esurn2(nolef1)+j;
                auxesurn1=esurn1(post_cont);
                M(lef, auxesurn1)=M(lef, auxesurn1)-mobility(ifacont)*normcont*parameter(1,1,ifacont)*w(post_cont);
            end
        end
        % contribuições do nó 2
        if nflag(nolef2,1)<200

            I(lef,1)=I(lef,1)+ mobility(ifacont)*normcont*parameter(1,2,ifacont)*nflag(nolef2,2);

        else

            for j=1:(esurn2(nolef2+1)-esurn2(nolef2))

                post_cont=esurn2(nolef2)+j;
                auxesurn1=esurn1(post_cont);
                M(lef, auxesurn1)=M(lef, auxesurn1)-mobility(ifacont)*normcont*parameter(1,2,ifacont)*w(post_cont);

            end
        end

    end
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxnolef1=parameter(1,3,ifactual);
    auxnolef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxnorel1=parameter(2,3,ifactual);
    auxnorel2=parameter(2,4,ifactual);
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*murel*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ARR=norma*mulef*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ mobility(ifactual)*ALL;
    M(lef,rel)=M(lef,rel)- mobility(ifactual)*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ mobility(ifactual)*ARR;
    M(rel,lef)=M(rel,lef)- mobility(ifactual)*ALL;
    % contribuições esquerda
    % termo 1
    if nflag(auxnolef1,1)>200
        
        for j=1:(esurn2(auxnolef1+1)-esurn2(auxnolef1))
            
            post_cont=esurn2(auxnolef1)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)-mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)+mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*w(post_cont);
        end
    else
        I(lef,1)=I(lef,1)+ mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
        I(rel,1)=I(rel,1)- mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnolef1,1)==202
        
        I(lef)=I(lef)+mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
        
        I(rel)=I(rel)-mobility(ifactual)*murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
    end
    % termo 2
    if nflag(auxnolef2,1)>200
        
        for j=1:(esurn2(auxnolef2+1)-esurn2(auxnolef2))
            
            post_cont=esurn2(auxnolef2)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)- mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)+ mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*w(post_cont);
        end
    else
        I(lef,1)=I(lef,1)+ mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        I(rel,1)=I(rel,1)- mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnolef2,1)==202
        
        I(lef)=I(lef)+mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
        
        I(rel)=I(rel)-mobility(ifactual)*murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
    end
    % contribuição do elemento a direita
    % termo 1
    if nflag(auxnorel1,1)>200
        
        for j=1:(esurn2(auxnorel1+1)-esurn2(auxnorel1))
            
            post_cont=esurn2(auxnorel1)+j;
            auxesurn1=esurn1(post_cont);
            M(lef, auxesurn1)=M(lef, auxesurn1)+mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*w(post_cont);
            M(rel, auxesurn1)=M(rel, auxesurn1)-mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*w(post_cont);

        end
    else
        I(lef,1)=I(lef,1)- mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
        I(rel,1)=I(rel,1)+ mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnorel1,1)==202
        
        I(lef)=I(lef)-mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
        
        I(rel)=I(rel)+mobility(ifactual)*mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
    end
    % termo 2
    if nflag(auxnorel2,1)>200
        
         for j=1:(esurn2(auxnorel2+1)-esurn2(auxnorel2))
             
             post_cont=esurn2(auxnorel2)+j;
             auxesurn1=esurn1(post_cont);
             M(lef, auxesurn1)=M(lef, auxesurn1)+ mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*w(post_cont);
             M(rel, auxesurn1)=M(rel, auxesurn1)- mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*w(post_cont);

         end
    else
        I(lef,1)=I(lef,1)- mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
        I(rel,1)=I(rel,1)+ mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
    end
    % no contorno de Neumann com fluxo diferente de zeros
    if nflag(auxnorel2,1)==202
        
        I(lef)=I(lef)-mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
        
        I(rel)=I(rel)+mobility(ifactual)*mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
    end
end

% adequação da matriz nos poços produtores
for iw=1:size(wells,1)      
    if (wells(iw,1)~=0)&&(wells(iw,5)>400)&&(wells(iw,5)<600) % Caso seja prescrita alguma pressão.
        M(wells(iw,1),:)=0*M(wells(iw,1),:);
        M(wells(iw,1),wells(iw,1))=1;
        I(wells(iw,1))=wells(iw,6);
    end
end

end

