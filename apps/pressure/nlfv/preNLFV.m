function [parameter,p_old,tol,nit,nflagface,auxface,gamma] = preNLFV(kmap,nflag)

global elem bcflag formethod bedge

auxface = 0;
nflagface = 0;

if strcmp(formethod,'MPFAQL')||strcmp(formethod,'NLFVPP')
    [parameter]=coefficientLPSangle(kmap,nflag);
else
    [faces] = element_face;
    % calculoa dos pontos armonicos
    [pointarmonic] = harmonicopoint(kmap);
    %calculo dos parametros ou constantes (ksi)
    [parameter,auxface] = coefficientPPSharmonicpoint(faces,pointarmonic,kmap);
    % adequa��o dos flag de face de contorno
    [nflagface]= contflagface(bedge);
end

g = min(bcflag(:,2));

% dados inicializa��o m�todos dos volumes finitos n�o linear
p_old=g*ones(size(elem,1),1);    % inicializando a press�o
tol=1e-8;                        % tolerancia
nit=1000;                        % # itera��es de Picard
gamma=0;

end