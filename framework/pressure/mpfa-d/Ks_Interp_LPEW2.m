function [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No, mobility,S,V,nw,no )
%Retorna os K(n ou t) necessários para a obtenção dos weights. kmap é a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);
global esurn2 esurn1 satlimit visc region 

nec=esurn2(No+1)-esurn2(No);

Kt1=zeros(nec,2); %As colunas representam i=1 e i=2.
Kt2=zeros(nec,1);
Kn1=zeros(nec,2);
Kn2=zeros(nec,1);
K=zeros(3);
K1=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

%Construção do tensor permeabilidade.%

%Cálculo das primeiras constantes, para todas as células que concorrem num%
%nó "ni".                                                                 %
for k=1:nec
        
j=esurn1(esurn2(No)+k);
   
    for i=1:2
        if (size(T,1)==size(O,1))&&(k==nec)&&(i==2)
            K(1,1)=mobility(V(i,k,No))*kmap(region(j),2);
            K(1,2)=mobility(V(i,k,No))*kmap(region(j),3);
            K(2,1)=mobility(V(i,k,No))*kmap(region(j),4);
            K(2,2)=mobility(V(i,k,No))*kmap(region(j),5);
            
            vt=T(1,:)-Qo;
            
            Kn1(k,i)=((R*(vt'))'*K*(R*(vt')))/norm(vt)^2;
            Kt1(k,i)=((R*(vt)')'*K*(vt)')/norm(vt)^2;
        else
            K(1,1)=mobility(V(i,k,No))*kmap(region(j),2);
            K(1,2)=mobility(V(i,k,No))*kmap(region(j),3);
            K(2,1)=mobility(V(i,k,No))*kmap(region(j),4);
            K(2,2)=mobility(V(i,k,No))*kmap(region(j),5);
            
             vt=T(k+i-1,:)-Qo;
            
            Kn1(k,i)=((R*(vt'))'*K*(R*(vt')))/norm(vt)^2;
            Kt1(k,i)=((R*(vt)')'*K*(vt)')/norm(vt)^2;
        end
    end
 
    %--------------- Calculo da mobilidade total  no elemento ------------%
    
    Krw1 = ((S(j) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;
    
    Kro1 = ((1 - S(j) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;
    
    L22 = Krw1/visc(1) +  Kro1/visc(2);   
    
    %------------------------- Tensores ----------------------------------%
    
    K1(1,1)= L22*kmap(region(j),2);
    K1(1,2)= L22*kmap(region(j),3);
    K1(2,1)= L22*kmap(region(j),4);
    K1(2,2)= L22*kmap(region(j),5);
    
    if (size(T,1)==size(O,1))&&(k==nec)
        
        %------------ Calculo dos K's internos no elemento ---------------%
        nvt=norm(T(1,:)-T(k,:));
        vt=T(1,:)-T(k,:);      
        
        Kn2(k)=((R*vt')'*K1*(R*vt'))/(nvt^2);
        Kt2(k)=((R*vt')'*K1*vt')/(nvt^2);
        
    else
        
        nvt=norm(T(k+1,:)-T(k,:));
        vt=T(k+1,:)-T(k,:);
        
        Kn2(k)=(R*vt')'*K1*(R*vt')/(nvt^2);
        Kt2(k)=((R*vt')'*K1*vt')/(nvt^2);
        
    end

end

end

