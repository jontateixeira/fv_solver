function[N,F,V,weightLS,esuel,esuel_coord,A,bound,S_old,S_cont]=presaturation(wells,nflag)

global satlimit elem centelem inedge
% adequação das faces por elemento e por nó em sentido anti-horario
% F: vetor de faces em cada elemento
% V: faces ordenados "convenientemente" ao rededor de um nó
% N: faces ordenados ao rededor de um no em sentido anti-horario
[F,V,N]=elementface(nflag,inedge);

[esuel, esuel_coord,A,bound] = esurnelem(inedge);

% calculo dos pesos para o método minimos quadrados
[weightLS] = weights(A,esuel,esuel_coord,centelem);

% Condiçao inicial saturaçao
S_old = zeros(size(elem,1),1);

S_old(:)=satlimit(2);
% Condiçao de contorno da saturaçao

S_cont=1-satlimit(1);

if isempty(wells)==1
    wells=0;
end

% adequação dos poços
if max(wells(:,1))~=0
    for i=1:size(wells,1)
        if wells(i,3)>300 && wells(i,3)<400
            S_old(wells(i,1))=S_cont;
        end
    end 
end

end