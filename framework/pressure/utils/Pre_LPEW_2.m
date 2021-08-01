function [w,s] = Pre_LPEW_2( kmap,mobility,V,S_old,nw,no,N )

global coord nsurn1 nsurn2 bcflag bedge inedge O1 P1 T1 ve21 ve11 theta21 theta11 neta1

% Retorna todos os parâmetros necessários às expressões dos fluxos.

apw=ones(size(coord,1),1);
r=zeros(size(coord,1),2);

for y=1:size(coord,1),
    
    No=y;
     
    Qo = coord(No,:);
    
    nq = sum(O1(:,1,y)~=-1);
    np = sum(P1(:,1,y)~=-1);
    
    O = O1(1:nq,1:3,y);
    P = P1(1:np,1:3,y);
    T = T1(1:np,1:3,y);
    ve2 = ve21(1:nq,y)+1e-12;
    ve1 = ve11(1:nq,y)+1e-12;
    theta2 = theta21(1:nq,y)+1e-12;
    theta1 = theta11(1:nq,y)+1e-12;
    neta = neta1(1:nq,1:2,y);
    
    % calculo dos Ks
    [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No,...
        mobility,S_old,V,nw,no );
       
    % calculo dos lamdas
    [ lambda,r ] = Lamdas_Weights_LPEW2( Kt1, Kt2, Kn1, Kn2, theta1,...
        theta2, ve1, ve2, neta, P, O,Qo,No,T,r);
    % calculo dos pesos
    for k=0:size(O,1)-1,
        w(apw(No)+k)=lambda(k+1)/sum(lambda); %Os pesos fazem sentido%%%%%%%%%%%
    end
    apw(No+1)=apw(No)+size(O,1);
    % interpolaçao das pressões nos contornos de Neumann
    vetor=nsurn1(nsurn2(No)+1:nsurn2(No+1));
    comp1=N(No,1);
    comp2=N(No,length(vetor));
    if comp1>size(inedge,1) && comp2>size(inedge,1)
        a=bcflag(:,1)==bedge(comp1-size(inedge,1),5);
        s1=find(a==1);
        b=bcflag(:,1)==bedge(comp2-size(inedge,1),5);
        s2=find(b==1);
        s(No,1)=-(1/sum(lambda))*(r(No,1)*bcflag(s1,2)+r(No,2)*bcflag(s2,2));
    end
    
end


end

