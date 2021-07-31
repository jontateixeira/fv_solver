function [ O, P, T, Qo ] = OPT_Interp_LPEW (No)

global coord esurn1 esurn2 nsurn1 nsurn2 centelem 
%Retorna os vetores O, P, T e Qo.
% Lembrando que estes esurn1, nsurn1 já estan ordenados em sentido
% anti-horario, sequencialmente. 

%Pré-alocação dos vetores.%

Qo=coord(No,:);                     % coordenada do nó "ni".

%Construção do vetor O, dos centróides (pontos de colocação) dos elementos
%que concorrem no nó ni.                                                  
O=centelem(esurn1(esurn2(No)+1:esurn2(No+1)),:);

%Construção dos vetores P, dos nós vizinhos ao nó "ni", e T, dos pontos
%médios das fases que concorrem no nó "ni".
P=coord(nsurn1(nsurn2(No)+1:nsurn2(No+1)),:);
T=0.5*(P+(ones(size(P,1),1)*Qo));

end
