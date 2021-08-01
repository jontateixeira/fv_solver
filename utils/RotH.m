function [RH]=RotH(vi)
%Função que retorna a rotação de um certo vetor em 90 graus, 
%anti-horariamente, na primeira coluna da matriz R, e horariamente, 
%na segunda. Restringe um vetor de 3 coordenadas a um de 2, considerando
%que a terceira coordenada é nula.
% vi2=zeros(2,1);
% if size(vi)~=[3 1]
%     vi=vi';
%     vi2(1)=vi(1);
%     vi2(2)=vi(2);
% end
RH=[0 1 0;-1 0 0;0 0 0]*vi';
end

