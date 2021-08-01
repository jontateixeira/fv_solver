function [mobility] = mobilityface(S_old,nw,no,S_cont)

global visc satlimit inedge bedge phflw

if phflw==2
    
    mobility=zeros(size(inedge,1)+size(bedge,1),1);

    for i = 1:size(inedge,1) % loop de fases internas para calcular as mobilidade e fluxo fracional

        ee=inedge(i,3);
        dd=inedge(i,4);

        Krw1 = ((S_old(ee) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;

        Kro1 = ((1 - S_old(ee) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;

        L1 = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo

        Krw2 = ((S_old(dd) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;

        Kro2 = ((1 - S_old(dd) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;

        L2 = Krw2/visc(1) +  Kro2/visc(2); % mobilidade total no elemento direito

        mobility(i)=(L1+L2)/2;

    end


    for i = 1:size(bedge,1) % loop de fases externas para calcular as mobilidade e fluxo fracional

        if bedge(i,5)==201 % 201 é o flag usado para fluxo igual a ZERO.

            ee=bedge(i,3);

            Krw1 = ((S_old(ee) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;

            Kro1 = ((1 - S_old(ee) - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;

            mobility(i+size(inedge,1)) = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo;

        else

            Krw1 = ((S_cont - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^nw;

            Kro1 = ((1 - S_cont - satlimit(1))/(1 - satlimit(1) - satlimit(2)))^no;

            mobility(i+size(inedge,1)) = Krw1/visc(1) +  Kro1/visc(2); % mobilidade total no elemento esquerdo

        end

    end

elseif phflw==1 || phflw==0
    
    mobility=ones(size(inedge,1)+size(bedge,1),1);
    
end


end
