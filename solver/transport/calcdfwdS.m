function [dfwdS] = calcdfwdS(Sw,nw,no)
%Define global parameters
global satlimit visc

%Define "dfdS" and "dgamadS" according to "dertype"
%Initialize some propoerties (two-phase flow)
Swi = satlimit(1);
Sor = satlimit(2);
miw = visc(1);
mio = visc(2);
%Initialize "dfdS" and "dgamadS"
dfwdS = zeros(length(Sw),1);
% dfodS = zeros(length(Sw),1);

%Fit parameter (water and oil)
kwmax = 1;
komax = 1;

%Calculate the derivate
for i = 1:length(Sw)
    %Define some terms:
    term1 = 1 - Swi - Sor;
    term2 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
    term3 = komax*((1 - ((Sw(i) - Swi)/term1))^no)/mio;
    term4 = kwmax*(((Sw(i) - Swi)/term1)^nw)/miw;
    term5 = kwmax*(((Sw(i) - Swi)/term1)^nw);
    term6 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1))/(term1*miw);
    term7 = no*komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1))/...
        (term1*mio);
    %Calculate "dfdS"
    dfwdS(i) = (term2/(term1*(term3 + term4)*miw)) - ...
        (term5*(term6 - term7))/(((term3 + term4)^2)*miw);
    
%     % Calculate dfodSw
%     
%     a = 1 - Swi - Sor;
%     u = (a^2 - 2*Sw(i)*Swi - (2*Sw(i) - 2*Swi)*(a^2) + Sw(i)^2 + Swi^2)/(mio*(a^2));
%     b = Sw(i)^2-2*Sw(i)*Swi+Swi^2;
%     c = Sw(i)^2+Swi^2-2*Sw(i)*Swi+2*Sw(i)-2*Swi;
%     
%     v = (b/(miw*a^2))+(1/mio)-(c/(mio*a^2));
%     dudSw = -(2*Swi - 2*Sw(i) + 2*(a^2))/(mio*(a^2));
%     dvdSw = (2/(a^2))*(((Sw(i)-Swi)/miw)-((Sw(i)-Swi+1)/mio));
%     
%     dfodS(i) = (v*dudSw - u*dvdSw)/(v^2); %regra do quociente
        
end

% dfdS = [dfwdS, dfodS];
% u = ((1-swi-sor)^2-[2*(sw-swi)*((1-swi-sor)^2)]+(sw^2-2*sw*swi+swi^2))/(((1-swi-sor)^2)*mio)
