function [Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap)

global bedge inedge coord centelem region
% Retorna alguns par�metros da express�o dos fluxos na face interna e no contorno.
% podemos verificar na pag. 5 e na pag. 6 do paper chin�s

%Prealoca��o das matrizes.%

Hesq=zeros(size(bedge,1),1);
Kde=zeros(size(inedge,1),1);
Ded=zeros(size(inedge,1),1);
Kn=zeros(size(bedge,1),1);
Kt=zeros(size(bedge,1),1);
% determinar o flag do n� interior e fronteira de Neumann
K1=zeros(3,3);
K2=zeros(3,3);
K=zeros(3,3);

%Loop de arestas de contorno.%

for ifacont=1:size(bedge,1)
    
    %Determina��o do baricentro a elemento � esquerda
    
    C1=centelem(bedge(ifacont,3),:);
    
    %Determina��o das alturas dos elementos � esquerda.
    
    ve1=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); % face
    ve2=coord(bedge(ifacont,2),:)-C1; %Do centro esquerdo ao fim da face.
    
    ce=cross(ve1,ve2); % produto vetorial
    Hesq(ifacont)=norm(ce)/norm(ve1); % altura a relativo as faces do contorno
    
    %Essa � UMA maneira de construir os tensores
    K(1,1)=kmap(region(bedge(ifacont,3)),2);
    K(1,2)=kmap(region(bedge(ifacont,3)),3);
    K(2,1)=kmap(region(bedge(ifacont,3)),4);
    K(2,2)=kmap(region(bedge(ifacont,3)),5);
    
    %C�lculo das constantes tangenciais e normais
    Kn(ifacont)=(RotH(ve1)'*K*RotH(ve1))/norm(ve1)^2;
    Kt(ifacont)=(RotH(ve1)'*K*(ve1)')/norm(ve1)^2;
    
    %Adequa��o dos flags dos n�s de Dirichlet
        
end
%Loop de arestas internas.%

for iface=1:size(inedge,1),
    
    %Determina��o dos centr�ides dos elementos � direita e � esquerda.%
    C1=centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
    C2=centelem(inedge(iface,4),:); % baricentro do elemento direito
    vcen=C2-C1;
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:); % ARESTA!
    nvd1=norm(vd1);
    
    %Determina��o das alturas dos centr�ides dos elementos � direita e �%
    %esquerda. 
    
    vd2=C2-coord(inedge(iface,1),:); %Do in�cio da aresta at� o centro da c�lula da direita.
    cd=cross(vd1,vd2);
    H2=norm(cd)/nvd1; % altura a direita
      
    ve2=C1-coord(inedge(iface,1),:);
    ce=cross(vd1,ve2);
    H1=norm(ce)/nvd1; % altura a esquerda
    
    %C�lculo das constantes.%
    %A segunda entrada ser� tal que: 1=dir, 2=esq.
    
    %Essa � UMA maneira de construir os tensores.
    
    K1(1,1)=kmap(region(inedge(iface,3)),2);
    K1(1,2)=kmap(region(inedge(iface,3)),3);
    K1(2,1)=kmap(region(inedge(iface,3)),4);
    K1(2,2)=kmap(region(inedge(iface,3)),5);
    
    % tensor a elemento a direita
    
    K2(1,1)=kmap(region(inedge(iface,4)),2);
    K2(1,2)=kmap(region(inedge(iface,4)),3);
    K2(2,1)=kmap(region(inedge(iface,4)),4);
    K2(2,2)=kmap(region(inedge(iface,4)),5);
    
    % calculo das constantes tangenciais e normais em cada face interna
    Kn1=(RotH(vd1)'*K1*RotH(vd1))/nvd1^2;
    Kt1=(RotH(vd1)'*K1*(vd1)')/nvd1^2;
    
    Kn2=(RotH(vd1)'*K2*RotH(vd1))/nvd1^2;
    Kt2=(RotH(vd1)'*K2*(vd1)')/nvd1^2;
    
    % calculo das constantes nas faces internas
    Kde(iface)=-nvd1*((Kn1*Kn2))/(Kn1*H2+Kn2*H1);
    
    % Ded: � uma constante que tem constantes geometricas + contantes
    % tangeciais
    Ded(iface)=(dot(vd1,vcen)/nvd1^2)-(1/nvd1)*((Kt2/Kn2)*H1+(Kt1/Kn1)*H2);
    
end
%------------------ Fim do loop de arestas internas ----------------------%


end

