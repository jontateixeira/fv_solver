function [faces, neigh] = teste_unified_preprocess

global elem inedge bedge 

nel = size(elem,1); 
nfel = 4;  
nedint = size(inedge,1);  % faces internas

%% Definindo faces de cada elemento ordenadas (w,e,s,n)
edges = [inedge(:,[1,2]); bedge(:,[1,2])];
ied = (1:size(edges))';

faces = zeros(nel, nfel); % faces dos elementos
lei = [1,2; 2,3; 3,4; 4,1];  % local edges index
for e=1:nel
    for i=1:nfel
        idx = ismember(edges,elem(e, lei(i,:)),'rows'); 
        s = sum(idx);
        if s==0
            idx = ismember(edges,elem(e, circshift(lei(i,:),[0,-1])),'rows');
            faces(e,i) = ied(idx);
        else
            faces(e,i) = s*ied(idx);
        end
    end
end
faces = faces(:,[4,2,1,3]); % ordenar -> (w,e,s,n)


%% Definindo elementos vizinhos ordenados (w,e,s,n)
neigh = zeros(size(faces)); % elementos vizinhos 
edges = [inedge(:,[3,4]); [bedge(:,3), zeros(size(bedge,1),1)]];
for e=1:nel
    for i=1:nfel
        ied = faces(e,i);
        neigh(e,i) = (ied<=nedint)*edges(ied,~ismember(edges(ied,:),e));
    end
end

end