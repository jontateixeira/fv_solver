function [ O1, P1, T1, ve21, ve11, theta21, theta11, neta1 ] = geomparamLPEW2( coord, esurn2 )
%
m = 0;
for i=2:size(esurn2)
    if (esurn2(i)-esurn2(i-1)) > m
        m = esurn2(i)-esurn2(i-1); 
    end
end

O1 = -1*ones(m,3,size(coord,1));
P1 = -1*ones(m,3,size(coord,1));
T1 = -1*ones(m,3,size(coord,1));
ve21 = -1*ones(m,size(coord,1));
ve11 = -1*ones(m,size(coord,1));
theta21 = -1*ones(m,size(coord,1));
theta11 = -1*ones(m,size(coord,1));
neta1 = -1*ones(m,2,size(coord,1));

for y=1:size(coord,1),
    
    No=y;
    % calculos dos vetores O, P, T, Q
    [ O, P, T, Qo ] = OPT_Interp_LPEW (No);
    % calculo dos angulos
    [ ve2, ve1, theta2, theta1 ] = angulos_Interp_LPEW2( O, P, T, Qo, No );
    % calculo dos netas
    [ neta ] = netas_Interp_LPEW( O, P, T, Qo, No );
    
    O1(1:size(O,1),1:3,y) = O;
    P1(1:size(P,1),1:3,y) = P;
    T1(1:size(T,1),1:3,y) = T;
    ve21(1:size(ve2,2),y) = ve2;
    ve11(1:size(ve1,2),y) = ve1;
    theta21(1:size(theta2,2),y) = theta2;
    theta11(1:size(theta1,2),y) = theta1;
    neta1(1:size(neta,1),1:2,y) = neta;
    
end

end

