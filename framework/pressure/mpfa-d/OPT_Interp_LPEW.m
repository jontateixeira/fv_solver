function [ O, P, T, Qo ] = OPT_Interp_LPEW (No)

global coord esurn1 esurn2 nsurn1 nsurn2 centelem 
%Retorna os vetores O, P, T e Qo.
% Lembrando que estes esurn1, nsurn1 j� estan ordenados em sentido
% anti-horario, sequencialmente. 

%Pr�-aloca��o dos vetores.%

Qo=coord(No,:);                     % coordenada do n� "ni".

%Constru��o do vetor O, dos centr�ides (pontos de coloca��o) dos elementos
%que concorrem no n� ni.                                                  
O=centelem(esurn1(esurn2(No)+1:esurn2(No+1)),:);

%Constru��o dos vetores P, dos n�s vizinhos ao n� "ni", e T, dos pontos
%m�dios das fases que concorrem no n� "ni".
P=coord(nsurn1(nsurn2(No)+1:nsurn2(No+1)),:);
T=0.5*(P+(ones(size(P,1),1)*Qo));

end
