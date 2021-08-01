function [ M, I ] = globalmatrix( w, s, Kde, Ded, Kn, Kt, nflag, Hesq,...
                                 wells, mobility)

global coord elem esurn1 esurn2 bedge inedge centelem elemarea bcflag...
       source

% Inicialização das matrizes ---------------------------------------------%
M=sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I=sparse(size(elem,1),1);
%-------------------------------------------------------------------------%

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

for ib=1:size(bedge,1)
    v0=coord(bedge(ib,2),:)-coord(bedge(ib,1),:);
    v1=centelem(bedge(ib,3),:)-coord(bedge(ib,1),:);
    v2=centelem(bedge(ib,3),:)-coord(bedge(ib,2),:);
    normcont=norm(coord(bedge(ib,2),:)-coord(bedge(ib,1),:));
    
    if bedge(ib,5)<200 % Caso haja arestas de Dirichlet
        if bedge(ib,1)~=bedge(ib,2)
            c1=nflag(bedge(ib,1),2);
            c2=nflag(bedge(ib,2),2);

            A=-Kn(ib)/(Hesq(ib)*norm(v0));

            M(bedge(ib,3),bedge(ib,3))=M(bedge(ib,3),...
                bedge(ib,3))-mobility(ib+size(inedge,1))*A*(norm(v0)^2);

            I(bedge(ib,3))=I(bedge(ib,3))-mobility(ib+...
                size(inedge,1))*(dot(v2,-v0)*c2+dot ...
                (v1,v0)*c1)*A+mobility(ib+size(inedge,1))*(c2-c1)*Kt(ib);
        else
            c1=nflag(bedge(ib,1),2);
            M(bedge(ib,3),bedge(ib,3))=M(bedge(ib,3),...
                bedge(ib,3))+mobility(ib+size(inedge,1))*Kn(ib)/Hesq(ib);
            I(bedge(ib,3))=I(bedge(ib,3))+mobility(ib+...
                size(inedge,1))*c1*Kn(ib)/Hesq(ib);
        end
    else
        x=bcflag(:,1)==bedge(ib,5);
        r=find(x==1);
        I(bedge(ib,3))=I(bedge(ib,3))-normcont*bcflag(r,2);
    end
end

% contribuição nas faces internas
for iface=1:size(inedge,1)        
    
    %Contabiliza as contribuições do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
   
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3),inedge(iface,3))...
        - mobility(iface)*Kde(iface);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3),inedge(iface,4))...
        + mobility(iface)*Kde(iface);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4),inedge(iface,4))...
        - mobility(iface)*Kde(iface);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4),inedge(iface,3))...
        + mobility(iface)*Kde(iface);
    
    %Se os nós das arestas estiverem em fronteiras de Dirichlet, suas
    %contribuições serão contabilizadas logo abaixo.
    
    if nflag(inedge(iface,1),1)<200
        I(inedge(iface,3))=I(inedge(iface,3))-mobility(iface)*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
        I(inedge(iface,4))=I(inedge(iface,4))+mobility(iface)*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
    end
    if nflag(inedge(iface,2),1)<200
        I(inedge(iface,3))=I(inedge(iface,3))+mobility(iface)*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
        I(inedge(iface,4))=I(inedge(iface,4))-mobility(iface)*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
    end
    
    % quando o nó pertece ao contorno de Neumann com fluxo prescrito
    % diferente de zero --------------------------------------------------%
    if (nflag(inedge(iface,1),1)>201)&&(nflag(inedge(iface,1),1)<300)
        
        I(inedge(iface,3))=I(inedge(iface,3))-mobility(iface)*Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))+mobility(iface)*Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
    end
    if (nflag(inedge(iface,2),1)>201)&&(nflag(inedge(iface,2),1)<300)
        
        I(inedge(iface,3))=I(inedge(iface,3))+mobility(iface)*Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))-mobility(iface)*Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
    end
    %---------------------------------------------------------------------%
    
    %Contabilização das contribuições dos nós que não estão na
    %fronteiras de Dirichlet.
    
    if nflag(inedge(iface,1),1)>200
        for j=1:(esurn2(inedge(iface,1)+1)-esurn2(inedge(iface,1)))
            
            post_cont=esurn2(inedge(iface,1))+j;
            
            M(inedge(iface,3), esurn1(post_cont))=M(inedge(iface,3),...
                esurn1(post_cont)) + mobility(iface)*Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont))=M(inedge(iface,4),...
                esurn1(post_cont)) - mobility(iface)*Kde(iface)*Ded(iface)*w(post_cont);
            
        end
    end
    if nflag(inedge(iface,2),1)>200
        for j=1:(esurn2(inedge(iface,2)+1)-esurn2(inedge(iface,2))),
            
            post_cont=esurn2(inedge(iface,2))+j;
            
            M(inedge(iface,3), esurn1(post_cont))=M(inedge(iface,3),...
                esurn1(post_cont)) - mobility(iface)*Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont))=M(inedge(iface,4),...
                esurn1(post_cont)) + mobility(iface)*Kde(iface)*Ded(iface)*w(post_cont);
        end
    end
   
end

for iw=1:size(wells,1)      
    if (wells(iw,1)~=0)&&(wells(iw,5)>400)&&(wells(iw,5)<600) % Caso seja prescrita alguma pressão.
        M(wells(iw,1),:)=0*M(wells(iw,1),:);
        M(wells(iw,1),wells(iw,1))=1;
        I(wells(iw,1))=wells(iw,6);
    end
end

end

