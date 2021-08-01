function [ netas ] = netas_Interp_LPEW( O, P, T, Qo, No )
%Retorna os netas.
global esurn2

%Prealocação do vetore.%

netas=zeros(esurn2(No+1)-esurn2(No),2);

%Loop que percorre os elementos em torno do nó "ni".%

for k=1:size(netas,1),
      
    %Preenchimento da segunda coluna do vetor "netas".%
    
    if (k==size(netas,1))&&(size(P,1)==size(O,1))
        v1=O(k,:)-Qo;
        v2=P(1,:)-Qo;
        ce=cross(v1,v2);
        tce=norm(ce);
        h2=tce/norm(v2);
        netas(k,2)=norm(T(1,:)-Qo)/h2;
    else
        v1=O(k,:)-Qo;
        v2=P(k+1,:)-Qo;
        ce=cross(v1,v2);
        tce=norm(ce);
        h2=tce/norm(v2);
        netas(k,2)=norm(T(k+1,:)-Qo)/h2;
    end
    
    %Preenchimento da primeira coluna do vetor "netas".%
    
    v1=O(k,:)-Qo;
    v2=P(k,:)-Qo;
    ce=cross(v1,v2);
    tce=norm(ce);
    h1=tce/norm(v2);
    netas(k,1)=norm(T(k,:)-Qo)/h1;
    
end

end

