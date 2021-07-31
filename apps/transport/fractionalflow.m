function [f_elem] = fractionalflow(S_old,nw,no)
    global visc satlimit phflw

    f_elem = S_old;
    if phflw==2
        % calculo das permeabilidade usando modelo Brooks e Corey
        % calculando as permeabilidades relativas e a mobilidade do elemento
        Krw2 = ((S_old - satlimit(1))/(1-satlimit(1)-satlimit(2))).^nw;
        Kro2 = ((1 - S_old - satlimit(1))/(1-satlimit(1)-satlimit(2))).^no;
        
        L2 = Krw2/visc(1) + Kro2/visc(2);
        f_elem = (Krw2/visc(1))./L2; 
    end
end
