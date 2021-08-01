function [local]=mapeamento_piola(current_cell, x_glob, y_glob, elem, coord)
%
%%Set coords function on reference space
% 
%% EXEMPLO RESOLVIDO  
% x1 = 9; x2 = 11; x3 = 11; x4 = 7; 
% y1 = 4; y2 = 2; y3 = 6; y4 = 7; 
% x_glob = 10;
% y_glob = 4.5;
% Coords on R = P(0.67,0.5); 
% 
% [local]=mapeamento_piola(C, elem, coord) 

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% original version
%Determinar tolerância
%tol = 1e-10;
% modified JCT
tol = 100*eps;
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

%% Nodes em P
for i = 1: size(elem(current_cell,:),1)
    
    
    %% Nodes em P
    x1=coord(elem(current_cell(i),1),1); x2=coord(elem(current_cell(i),2),1); 
    x3=coord(elem(current_cell(i),3),1); x4=coord(elem(current_cell(i),4),1);
    y1=coord(elem(current_cell(i),1),2); y2=coord(elem(current_cell(i),2),2); 
    y3=coord(elem(current_cell(i),3),2); y4=coord(elem(current_cell(i),4),2);

    %% Inverse bilinear transform (P -> R)

    %Reference: 2003_Haegland_Streamline_Tracing_on_Irregular_Grids
    %Appendix_E(pag. 96)
    ax = x1 - x2 + x3 - x4;        
    bx = x1 - x2;
    cx = x1 - x4;
    dx = x_glob(i,1) - x1;
    ay = y1 - y2 + y3 - y4;
    by = y1 - y2;
    cy = y1 - y4;
    dy = y_glob(i,1) - y1;

    ex = x3 - x4;
    ey = y3 - y4; 

    r = cx*ay - cy*ax;
    s = bx*ay - by*ax; %s= bx*ey - by*ex; 
    t = dx*ay - dy*ax; 

    if ax ~=0 && s~=0 && r~=0 %(CASE 1)

        A = (-r*ax);
        B = (r*bx - (t*ax) -(cx*s));
        C = (t*bx - (dx*s));

        % Coord Y on R       
        y_ref(1) = (-(B) + (sqrt(B^2-4*A*C)))/(2*A);
        y_ref(2) = (-(B) - (sqrt(B^2-4*A*C)))/(2*A);
                     
        y_ref(abs(abs(y_ref)-1)<tol)=round(y_ref(abs(abs(y_ref)-1)<tol));
        
        y_ref(abs(y_ref)<tol)=round(y_ref(abs(y_ref)<tol));
        
        if y_ref(1)>=0 && y_ref(1)<=1
            y_local(i) = y_ref(1);
        elseif y_ref(2)>=0 && y_ref(2)<=1
            y_local(i) = y_ref(2);
        else
            y_aux = 10e5*ones(1,4);
            if y_ref(1)>1
                y_aux(1) = y_ref(1)-1;
            elseif y_ref(1)<0
                y_aux(2) = 0-y_ref(1);
            end
            if y_ref(2)>1
                y_aux(3) = y_ref(2)-1;
            elseif y_ref(2)<0
                y_aux(4) = 0-y_ref(2);
            end
            y_dis = find(y_aux == min(y_aux));
            if y_dis == 1 || y_dis == 2
                y_ref(1) = round (y_ref(1));
                y_local(i) = y_ref(1);
            elseif y_dis == 3 || y_dis == 4
                y_ref(2) = round (y_ref(2));
                y_local(i) = y_ref(2);
            end
        end
        

%         y_local(i) = y_ref((y_ref>=0)&(y_ref<=1));
        
        % Coord X on R
        x_local(i) = (-(r*y_local(i)/s) - (t/s));
       
    end
        
    if ax ~=0 && s~=0 && r == 0 %(CASE 2)
        x_local(i) = -(t/s);
        y_local(i) = -((t*bx - dx*s)/(t*ax +cx*s));
    end

    if ax ~=0 && s == 0 && r~=0 %(CASE 3)
        x_local(i) = -((cx*t - dx*r)/(r*bx +ax*t));
        y_local(i) = -(t/r);
    end

    if ax ==0 && cx ~=0 %(CASE 4)
        
        if ay ==0 || bx ==0
           
            x_local(i) = (dy*cx - cy*dx)/(cy*bx -by*cx -ay*dx);
            y_local(i) = ((-bx*x_local(i))-dx)/cx;
            
        elseif ay ~=tol && bx ~=tol
        
        A = (-bx*ay);
        B = (cx*bx - by*cx - ay*cx);
        C = (cy*dx - dy*cx);

        % Coord X on R     
        
        x_ref(1) = (-(B) + (sqrt(B^2-4*A*C)))/(2*A);
        x_ref(2) = (-(B) - (sqrt(B^2-4*A*C)))/(2*A);
        
        x_ref(abs(abs(x_ref)-1)<tol)=round(x_ref(abs(abs(x_ref)-1)<tol));
        
        x_ref(abs(x_ref)<tol)=round(x_ref(abs(x_ref)<tol));
                    
        if x_ref(1)>=0 && x_ref(1)<=1
            x_local(i) = x_ref(1);
        elseif x_ref(2)>=0 && x_ref(2)<=1
            x_local(i) = x_ref(2);
        else
            x_aux = 10e5*ones(1,4);
            if x_ref(1)>1
                x_aux(1) = x_ref(1)-1;
            elseif x_ref(1)<0
                x_aux(2) = 0-x_ref(1);
            end
            if x_ref(2)>1
                x_aux(3) = x_ref(2)-1;
            elseif x_ref(2)<0
                x_aux(4) = 0-x_ref(2);
            end
            x_dis = find(x_aux == min(x_aux));
            if x_dis == 1 || x_dis == 2
                x_ref(1) = round (x_ref(1));
                x_local(i) = x_ref(1);
            elseif x_dis == 3 || x_dis == 4
                x_ref(2) = round (x_ref(2));
                x_local(i) = x_ref(2);
            end
        end
             
%         x_local(i) = x_ref((x_ref>=0)&(x_ref<=1));

        % Coord Y on R
        y_local(i) = ((-bx*x_local(i))-dx)/cx;
        
        end
        
    end

    if ax ==0 && cx ==0 %(CASE 5)
        x_local(i) = -(dx/bx);
        y_local(i) = (by*dx - dy*bx)/(ay*dx +cy*bx);
    end


end


y_local(abs(abs(y_local)-1)<tol)=round(y_local(abs(abs(y_local)-1)<tol));
y_local(abs(y_local)<tol)=round(y_local(abs(y_local)<tol));

x_local(abs(abs(x_local)-1)<tol)=round(x_local(abs(abs(x_local)-1)<tol));
x_local(abs(x_local)<tol)=round(x_local(abs(x_local)<tol));

 local = [current_cell, x_local', y_local'];

%% TESTE

 local(local(:,2)<tol,2) = 0; local(local(:,end)<tol,end) = 0;
 local(local(:,2)>.999999,2) = 1; local(local(:,end)>0.999999,end) = 1;


end
