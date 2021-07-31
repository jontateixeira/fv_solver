function [M,I]=assemblematrixDMPSY(p,pinterp,gamma,parameter,weightDMP,auxface,wells,mobility)

global inedge coord bedge bcflag elem elemarea source

valuemin=1e-16;
I=sparse(size(elem,1),1);
M=sparse(size(elem,1),size(elem,1));
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;
%% fonte
I=I+source.*elemarea;

if max(wells)~=0
    sumvol=0;
    for iw = 1:size(wells,1)
        
        if wells(iw,2)==1            % injetor
            I(wells(iw,1))= 1*elemarea(wells(iw,1));        % injeta um m3 de agua por dia (d)
            sumvol=sumvol+ elemarea(wells(iw,1));
        end
    end
    I=I./sumvol;
else
    for ifacont=1:size(bedge,1)
        lef=bedge(ifacont,3);
        
        normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
        
        if bedge(ifacont,5)>200
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            I(lef)=I(lef)- normcont*bcflag(r,2);
        else
            % contriuição no elemento do contorno
            
            [Transmicont1,elemavaliar1,pressureface1,contorno1]=auxassemblematrixcontourDMPSY(parameter(1,3,ifacont),...
                parameter(1,1,ifacont),lef,pinterp,normcont,weightDMP,bedge);
            
            %auxelemavaliar1=auxface(ifacont,1);
            
            M(lef,elemavaliar1)=M(lef,elemavaliar1)-(1-contorno1)*Transmicont1*mobility(ifacont,1);
            
            [Transmicont2,elemavaliar2,pressureface2,contorno2]=auxassemblematrixcontourDMPSY(parameter(1,4,ifacont),...
                parameter(1,2,ifacont),lef,pinterp,normcont,weightDMP,bedge);
            
            %auxelemavaliar2=auxface(ifacont,2);
            
            M(lef,elemavaliar2)=M(lef,elemavaliar2)-(1-contorno2)*Transmicont2*mobility(ifacont,1);
            
            %__________________________________________________________________
            
            M(lef,lef)=M(lef,lef)+ Transmicont2 + Transmicont1*mobility(ifacont,1);
            
            I(lef)=I(lef) + contorno1*Transmicont1*pressureface1*mobility(ifacont,1) + contorno2*Transmicont2*pressureface2*mobility(ifacont,1);
            %__________________________________________________________________
        end
        
    end
end
x=1000;
y=1000;
% Montagem da matriz global

for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % calculo da norma do vetor "vd1"
    norma=norm(vd1);
    % Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    
    %-----------------------------------------------------------------%
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a esquerda
    
    ifacelef1=parameter(1,3,ifactual);
    ifacelef2=parameter(1,4,ifactual);
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a direita
    ifacerel1= parameter(2,3,ifactual);
    ifacerel2= parameter(2,4,ifactual);
    % Fluxo parcial do elemento a esquerda
    
    F1=norma*(parameter(1,1,ifactual)*(p(lef)-pinterp(ifacelef1))+ parameter(1,2,ifactual)*(p(lef)-pinterp(ifacelef2)));
    % Fluxo parcial do elemento a direita
    F2=norma*(parameter(2,1,ifactual)*(p(rel)-pinterp(ifacerel1))+ parameter(2,2,ifactual)*(p(rel)-pinterp(ifacerel2)));
    if abs(F1)<1e-20
        F1=0;
    end
    if abs(F2)<1e-20
        F2=0;
    end
    % caso 1 F1*F2<=0  pag. 2594 Sheng and yuan, 2011
    if F2*F1<0 || F2*F1==0
        % caso 1.1 e caso 1.2, idea adoptado de Gao and Wu, 2013
        mu1=(abs(F2)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        mu2=(abs(F1)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        betalef=mu1*(1-sign(F1*F2));
        betarel=mu2*(1-sign(F1*F2));
        
        % contribuições do elemento a esquerda
        % somando 1
        
        [flaglef11,flaglef12,auxweightlef1,pressureface1]=auxassemblematrixintDMPSY1(ifacelef1,lef,pinterp,weightDMP,bedge);
        % se contornolef1 é zero então o ifacelef1 pertece ao face interior
        % da malha
        
        auxTransmilef1= flaglef11*auxweightlef1*betalef*parameter(1,1,ifactual)*norma + flaglef12*betalef*parameter(1,1,ifactual)*norma;
        
        
        M(lef,auxface(ifactual,1))=  M(lef,auxface(ifactual,1))- flaglef11*auxTransmilef1*mobility(ifactual);
        
        % somando 2
        [flaglef21,flaglef22,auxweightlef2,pressureface2]=auxassemblematrixintDMPSY1(ifacelef2,lef,pinterp,weightDMP,bedge);
        
        auxTransmilef2= flaglef21*auxweightlef2*betalef*parameter(1,2,ifactual)*norma+ flaglef22*betalef*parameter(1,2,ifactual)*norma;
        
        
        M(lef,auxface(ifactual,2))=  M(lef,auxface(ifactual,2))- flaglef21*auxTransmilef2*mobility(ifactual);
        
        %__________________________________________________________________
        
        M(lef,lef)=M(lef,lef) + auxTransmilef2 + auxTransmilef1*mobility(ifactual);
        
        I(lef)=I(lef)+ flaglef22*auxTransmilef2*pressureface2*mobility(ifactual) + flaglef12*auxTransmilef1*pressureface1*mobility(ifactual);
        %__________________________________________________________________
        
        % contribuições do elemento a direita
        % somando 1
        
        [flagrel11,flagrel12,auxweightrel1,pressurefacerel1]=auxassemblematrixintDMPSY1(ifacerel1,rel,pinterp,weightDMP,bedge);
        
        auxTransmirel1= flagrel11*auxweightrel1*betarel*parameter(2,1,ifactual)*norma + flagrel12*betarel*parameter(2,1,ifactual)*norma;
        
        
        M(rel,auxface(ifactual,3))=  M(rel,auxface(ifactual,3))- flagrel11*auxTransmirel1*mobility(ifactual);
        
        % somando 2
        
        [flagrel21,flagrel22,auxweightrel2,pressurefacerel2]=auxassemblematrixintDMPSY1(ifacerel2,rel,pinterp,weightDMP,bedge);
        
        auxTransmirel2= flagrel21*auxweightrel2*betarel*parameter(2,2,ifactual)*norma + flagrel22*betarel*parameter(2,2,ifactual)*norma;
        
        
        M(rel,auxface(ifactual,4))=  M(rel,auxface(ifactual,4))- flagrel21*auxTransmirel2*mobility(ifactual);
        
        %__________________________________________________________________
        
        M(rel,rel)=M(rel,rel)+ auxTransmirel1*mobility(ifactual) + auxTransmirel2*mobility(ifactual);
        
        I(rel)=I(rel) + flagrel11*auxTransmirel1*pressurefacerel1*mobility(ifactual) + flagrel22*auxTransmirel2*pressurefacerel2*mobility(ifactual);
        %__________________________________________________________________
        
        % caso 2 F1*F2<=0  pag. 2594 Sheng and yuan, 2011
    else
        [F1b,]=calfluxopartial2(ifacelef1,ifacelef2,parameter(1,1,ifactual), parameter(1,2,ifactual), gamma,lef,pinterp,ifactual,p,norma);
        % aveliando no elemento a direita
        
        [F2b,]=calfluxopartial2(ifacerel1,ifacerel2,parameter(2,1,ifactual), parameter(2,2,ifactual), gamma,rel,pinterp,ifactual,p,norma);
        
        
        % faço esta escolha ja que um das faces pode coincidir com a face
        % em questão
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
            
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
            
        else
            % pode ser o caso que a face em questão seja nenhum dos outros
            % faces de interpolação
            auxparameter1=0;
            y=0;
            
        end
        %-----------------------------------------------------------------%
        
        % faço esta escolha ja que um das faces pode coincidir com a face
        % em questão
        if ifacerel1==ifactual
            auxparameter2=parameter(2,1,ifactual);
            
        elseif ifacerel2==ifactual
            
            auxparameter2=parameter(2,2,ifactual);
            
        else
            % pode ser o caso que a face em questão seja nenhum dos outros
            % faces de interpolação
            auxparameter2=0;
            x=0;
            
        end
        
        % Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
        % veja a REMARK 3.2 pag. 313 do mesmo artigo
        if F1b*F2b>0
            mulef=abs(F2b)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=abs(F1b)/(abs(F1b)+abs(F2b)+2*valuemin);
        else
            mulef=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        end
        % Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
        
        mulef1=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        murel1=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        % implementação do beta
        betalef=mulef*(1-sign(F1b*F2b));
        betarel=murel*(1-sign(F1b*F2b));
        %=================================================================%
        % caso 2.1
        if F1b*F2b<0 || F1b*F2b==0 || (F2b<0 && F1b<0 && x==0 && y==0) || (F2b>0 && F1b>0 && x==0 && y==0)
            % Calculo das contribuições do elemento a esquerda
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            
            % contribuição da transmisibilidade no elemento esquerda
            
            % somando 1
            
            [flaglef11,flaglef12,flaglef13,auxweightlef1,pressurefacelef1]=auxassemblematrixint(ifacelef1,lef,pinterp,ifactual,weightDMP,bedge);
            
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacelef1==ifactual
            Transmilef1=flaglef11*weightrel*gamma*betalef*parameter(1,1,ifactual)*norma;
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacelef1~=ifactual
            auxTransmilef1= flaglef12*auxweightlef1*betalef*parameter(1,1,ifactual)*norma + flaglef13*betalef*parameter(1,1,ifactual)*norma;
            
            M(lef,auxface(ifactual,1))=M(lef,auxface(ifactual,1))- flaglef12*auxTransmilef1*mobility(ifactual);
            
            % somando 2
            
            [flaglef21,flaglef22,flaglef23,auxweightlef2,pressurefacelef2]=auxassemblematrixint(ifacelef2,lef,pinterp,ifactual,weightDMP,bedge);
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacelef2==ifactual
            Transmilef2=flaglef21*weightrel*gamma*betalef*parameter(1,2,ifactual)*norma;
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacelef2~=ifactual ou quando a face
            % pertenece ao contorno.
            auxTransmilef2= flaglef22*auxweightlef2*betalef*parameter(1,2,ifactual)*norma  + flaglef23*betalef*parameter(1,2,ifactual)*norma;
            
            M(lef,auxface(ifactual,2))=M(lef,auxface(ifactual,2))-flaglef22*auxTransmilef2*mobility(ifactual);
            
            %______________________________________________________________
            
            M(lef,lef)=M(lef,lef) + mobility(ifactual)*(Transmilef2 + auxTransmilef2 + Transmilef1+ auxTransmilef1 + alfa);
            
            M(lef,rel)=M(lef,rel) - Transmilef1*mobility(ifactual) - Transmilef2*mobility(ifactual) - alfa*mobility(ifactual);
            
            I(lef)=I(lef) + flaglef13*auxTransmilef1*pressurefacelef1*mobility(ifactual) + flaglef23*auxTransmilef2*pressurefacelef2*mobility(ifactual);
            %______________________________________________________________
            
            % contribuição da transmisibilidade no elemento direita
            
            % somando 1
            [flagrel11,flagrel12,flagrel13,auxweightrel1,pressurefacerel1]=auxassemblematrixint(ifacerel1,rel,pinterp,ifactual,weightDMP,bedge);
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacerel1==ifactual
            Transmirel1=flagrel11*weightlef*gamma*betarel*parameter(2,1,ifactual)*norma;
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacerel1~=ifactual
            auxTransmirel1= flagrel12*auxweightrel1*betarel*parameter(2,1,ifactual)*norma + flagrel13*betarel*parameter(2,1,ifactual)*norma;
            
            M(rel,auxface(ifactual,3))=M(rel,auxface(ifactual,3))-flagrel12*auxTransmirel1*mobility(ifactual);
            
            % somando 2
            
            [flagrel21,flagrel22,flagrel23,auxweightrel2,pressurefacerel2]=auxassemblematrixint(ifacerel2,rel,pinterp,ifactual,weightDMP,bedge);
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacerel2==ifactual
            Transmirel2=flagrel21*weightlef*gamma*betarel*parameter(2,2,ifactual)*norma;
            % esta transmisibilidade é ativado quando ifacelef1 pertenece
            % ao interior da malha e ifacerel2~=ifactual ou quando a face
            % pertenece ao contorno.
            auxTransmirel2= flagrel22*auxweightrel2*betarel*parameter(2,2,ifactual)*norma + flagrel23*betarel*parameter(2,2,ifactual)*norma;
            
            M(rel,auxface(ifactual,4))=M(rel,auxface(ifactual,4))-flagrel22*auxTransmirel2*mobility(ifactual);
            
            %______________________________________________________________
            
            M(rel,rel)= M(rel,rel) + mobility(ifactual)*(Transmirel1 + auxTransmirel1 + Transmirel2 + auxTransmirel2 + alfa);
            
            M(rel,lef)= M(rel,lef) - Transmirel1*mobility(ifactual) - Transmirel2*mobility(ifactual) - alfa*mobility(ifactual);
            
            I(rel)=I(rel) + flagrel13*auxTransmirel1*pressurefacerel1*mobility(ifactual) + flagrel23*auxTransmirel2*pressurefacerel2*mobility(ifactual);
            
            %=============================================================%
            %  caso 2.2
        else
            
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            % Calculo da transmisibilidade (eq. 18)
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            % contribuição da transmisibilidade no elemento esquerda
            M(lef,lef)=M(lef,lef)+ alfa*mobility(ifactual);
            M(lef,rel)=M(lef,rel)- alfa*mobility(ifactual);
            % contribuição da transmisibilidade no elemento direita
            M(rel,rel)=M(rel,rel)+ alfa*mobility(ifactual);
            M(rel,lef)=M(rel,lef)- alfa*mobility(ifactual);
        end
        
    end
    
end
if max(wells)~=0
    for iw = 1:size(wells,1)
        if wells(iw,2)==2 %produtor
            M(wells(iw,1),:)=0*M(wells(iw,1),:);
            M(wells(iw,1),wells(iw,1))=1;
            I(wells(iw,1))=0;
        end
    end
end
% switch benchmark
%     case 'gaowu6'
%         % malha 23x23
% %          M(357,:)=0*M(357,:);
% %          M(357,357)=1;
% %          I(357)=1;
% %          M(173,:)=0*M(173,:);
% %          M(173,173)=1;
% %          I(173)=0;
%         % malha 11x11
%         %M(83,:)=0*M(83,:);
%         %M(83,83)=1;
%         %I(83)=1;
%         %M(39,:)=0*M(39,:);
%         %M(39,39)=1;
%         %I(39)=0;
% end
end