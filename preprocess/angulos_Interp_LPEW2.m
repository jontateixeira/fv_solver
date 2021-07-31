function [ ve2, ve1, theta2, theta1 ] = angulos_Interp_LPEW2( O, P, T, Qo,No )
%
global esurn2 

ve2=zeros(1,esurn2(No+1)-esurn2(No)); 
ve1=zeros(1,esurn2(No+1)-esurn2(No)); 
theta2=zeros(1,esurn2(No+1)-esurn2(No));
theta1=zeros(1,esurn2(No+1)-esurn2(No)); 

for k=1:size(ve2,2),
    
    %Determinação dos vetores necessários à obtenção dos cossenos:
    v0=O(k,:)-Qo;
    if (k==size(ve2,2))&&(size(P,1)==size(O,1))
        vth2=T(1,:)-Qo;
        v1=T(1,:)-T(k,:);
        nv1=norm(v1);
    else
        v1=T(k+1,:)-T(k,:);
        nv1=norm(v1);
        vth2=T(k+1,:)-Qo;
    end
    vth1=T(k,:)-Qo;
    %Determinação dos ângulos:
    ve1(k)=real(acos(dot(-vth1,v1)/(nv1*norm(vth1))));   % revisar esses signos
    ve2(k)=real(acos(dot(-vth2,-v1)/(nv1*norm(vth2))));
    theta2(k)=real(acos(dot(v0,vth2)/(norm(v0)*norm(vth2))));
    theta1(k)=real(acos(dot(v0,vth1)/(norm(v0)*norm(vth1))));

end

end

