function [ pT, pV ] = calcpT( p, mobility, S_old, V, nw, no, w )
% calculate the dynamic points values

global O1 P1 T1 ve21 ve11 theta21 theta11 neta1 coord inedge kmap esurn2 ...
       esurn1

pV = zeros(size(coord,1),1);
pT = zeros(size(inedge,1),1);

for y=1:size(coord,1),
    
    No=y;
    
    Qo = coord(No,:);
    
    nq = sum(O1(:,1,y)~=-1);
    np = sum(P1(:,1,y)~=-1);
    
    O = O1(1:nq,1:3,y);
    P = P1(1:np,1:3,y);
    T = T1(1:np,1:3,y);
    ve2 = ve21(1:nq,y);
    ve1 = ve11(1:nq,y);
    theta2 = theta21(1:nq,y);
    theta1 = theta11(1:nq,y);
    neta = neta1(1:nq,1:2,y);
    
    % calculo dos Ks
    [ Kt1, Kt2, Kn1, Kn2 ] = Ks_Interp_LPEW2( O, T, Qo, kmap, No,...
        mobility,S_old,V,nw,no );
    
    nec=size(O,1);
    vel = esurn1(esurn2(No)+1:esurn2(No+1));
    vwe = w(esurn2(No)+1:esurn2(No+1));
    vwe = vwe';
    vpe = zeros(1,size(vwe,1));
    for i=1:nec
        vpe(i) = p(vel(i));
    end
    pV(y) = vpe*vwe;
       
    if nq==np
        a=1;
    else
        a=2;
    end
    
    for k=a:nec,
        % Se for usar a equação dos pesos --------------------------------%
        if k==1
            c=nec;
        else
            c=k-1;
        end
        eK = vel(k);
        e1K = vel(c);
        for j=1:size(inedge,1)
           if (inedge(j,3)==eK && inedge(j,4)==e1K)||(inedge(j,4)==eK && inedge(j,3)==e1K)
               f=j;
           end
        end
        pK = vpe(k);
        p1K = vpe(c);
        %-----------------------------------------------------------------%
        
        pT(f) = pV(y) + ((Kn1(k,1)*neta(k,1)*(pK-pV(y))+Kn1(c,2)*neta(c,2)*(p1K-pV(y)))/(Kn1(c,2)*cot(theta2(c))+Kn1(k,1)*cot(theta1(k))-Kt1(c,2)+Kt1(k,1)));
        
    end
        
end

end

