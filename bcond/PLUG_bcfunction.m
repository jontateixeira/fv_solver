%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 30/06/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the value to boundary condition according the 
%benchmark which intend to run. 

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

function [bcattrib] = PLUG_bcfunction(vertices,flagptr)
%Define global parameters
global coord bcflag

%Define a boolean parameter
boolean = length(vertices) == 1;
numcase = 99;
%Define the midpoint coordinate
coordmid = mean(coord(vertices,1:2))*(1 - boolean) + ...
    coord(vertices(1),1:2)*boolean;
%Calculate the value of the boundary condition to parameters attributed
%throughout the function
    switch numcase
        %------------------------------------------------------------------
        %Example 2: homogeneous media with diagonal permeability tensor and 
        %Dirichlet boundary condition. Edwards and M. Pal, 2008 (Case 2)  
        case 2
            bcattrib = cos(pi*coordmid(1))*cosh(pi*coordmid(2));

        %------------------------------------------------------------------
        %Example 3: This case is a quadratic variation of pressure in a two
        %material domain. The permeabilities are orthotropic and the
        %boundary condition applied are of Dirichlet type. (Edwards and 
        %Zheng, 2008)
        case 3
            %Definition of parameters:
            a = 1/50;
            b = 1/10;
            ar = 1;
            f = 4*ar/(((a - 2)*b) + 1);
            br = (b - 1)*f;
            cr = f;
            dr = -cr*1/10;
            cl = a*b*cr;
            dl = dr;
            
            %When the domain is minor than 1/2 in the direction x
            if coordmid(1) < 0.5
                bcattrib = (cl*(coordmid(1)^2)) + (dl*(coordmid(2)^2));
            %When the domain is minor than 1/2 in the direction x
            elseif coordmid(1) >= 0.5
                bcattrib = ar + (br*coordmid(1)) + (cr*(coordmid(1)^2)) + ...
                    (dr*(coordmid(2)^2));
            end  %End of IF
                
        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 2). 
        %Equation 4.4 (alfa, ai and bi variations) - Two materials with 
        %k1 = 100 and k2 = 1 (angle = 2*pi/3)
        case 4.1
            %Definition of parameters
            alfa = 0.75472745;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 1;
            a1 = 1; 
            a2 = 100.980198;
            a3 = 100.980198;
            a4 = 100.980198;
            b1 = 1.00995049;
            b2 = 1.99990197;
            b3 = 1.99990197;
            b4 = 1.99990197;
            %Calculate the pressure field
            bcattrib = ...
                calcCase4(coordmid,numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4);

        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 3). 
        %Equation 4.4 (alfa, ai and bi variations) - Four materials: 
        %k1 = k3 = 100; k2 = k4 = 1 (first example of apendix, angle 2*pi/3)
        case 4.2
            %Definition of parameters
            alfa = 0.13448835;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 1;
            a1 = 1;
            a2 = 4.90138222;
            a3 = -0.85392910;
            a4 = -9.94074425;
            b1 = 0.14177447;
            b2 = -13.3407815;
            b3 = -0.53935618;
            b4 = 10.1578346;
            %Calculate the pressure field
            bcattrib = ...
                calcCase4(coordmid,numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4);

        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 4). 
        %Equation 4.4 (alfa, ai and bi variations) - Four materials: 
        %k1 = k3 = 5; k2 = k4 = 1, angle = pi/2
        case 4.3
            %Definition of parameters
            alfa = 0.53544095;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 2;
            a1 = 1;
            a2 = 2.33333333;
            a3 = 0.55555556;
            a4 = -0.48148148;
            b1 = 0.44721360;
            b2 = -0.74535599;
            b3 = -0.94411759;
            b4 = -2.40170264;
            %Calculate the pressure field
            bcattrib = ...
                calcCase4(coordmid,numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4);

        %------------------------------------------------------------------
        %Zheng Thesis, pp. 75 - Four materials: 
        %k1 = k3 = 6; k2 = k4 = 1, angle = pi/3
        case 4.4
            %Definition of parameters
            alfa = 0.51671199;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 3;
            %Other parameters
            a1 = 0.27735010;
            a2 = -0.91129318;
            a3 = -0.98406726;
            a4 = -1.75974652;
            b1 = 1;
            b2 = 1.71428571;
            b3 = 0.32944606;
            b4 = -0.820074971;
            %Calculate the pressure field
            bcattrib = ...
                calcCase4(coordmid,numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4);
            
        %------------------------------------------------------------------
        %Example 5.1: Aavatsmark and Eigestad, 2006 (see eqaution 33).
        %Anysotrppi ratio: 3 (?)
        case 5.1
            %Definition of important parameters:
            kapla = 3;
            %Calculate the pressure fields
            bcattrib = calcCase5(coordmid,kapla);

        %------------------------------------------------------------------
        %Example 5.2: Aavatsmark and Eigestad, 2006 (see eqaution 33).
        %Anysotrppi ratio: 10 (?)
        case 5.2
            %Definition of important parameters:
            kapla = 10;
            %Calculate the pressure fields
            bcattrib = calcCase5(coordmid,kapla);

        %------------------------------------------------------------------
        %Example 8.3: Aavatsmark and Eigestad, 2006 (see eqaution 33).
        %Anysotrppi ratio: 100 (?)
        case 5.3
            %Definition of important parameters:
            kapla = 100;
            %Calculate the pressure fields
            bcattrib = calcCase5(coordmid,kapla);

        %------------------------------------------------------------------
        %Example 6.1: Aavatsmark and Eigestad, 2006 (see equation 32).
        %Anysotropi ratio 1/1000
        case 6.1
            %Definition of parameters
            kapla = 1e-3;
            %Calculate the presure over domain's boundary
            bcattrib = calcCase6(coordmid,kapla);

        %------------------------------------------------------------------
        %Example 6.2: Aavatsmark and Eigestad, 2006 (see equation 32).
        %Anysotropi ratio 100
        case 6.2
            %Definition of important parameters:
            kapla = 100;
            %Calculate the presure over domain's boundary
            bcattrib = calcCase6(coordmid,kapla);
         
        case 7.1
            x=coordmid(1);
            y=coordmid(2);
            
            if (((0<x || x==0 )&& (x <0.2 || x==0.2)) && y==0) || (((0<y || y==0 )&& (y <0.2 || y==0.2)) && x==0)
               
                bcattrib=1;
            elseif (((0.8<x || x==0.8 )&& (x <1 || x==1)) && y==1) || (((0.8<y || y==0.8 )&& (y <1 || y==1)) && x==1)
               
                bcattrib=0;
            elseif (((0.3<x || x==0.3 )&& (x <1 || x==1)) && y==0) || (((0.3<y || y==0.3 )&& (y <1 || y==1)) && x==0)
                
                bcattrib=0.5;
            elseif (((0<x || x==0 )&& (x <0.7 || x==0.7)) && y==1) || (((0<y || y==0 )&& (y <0.7 || y==0.7)) && x==1)
                
                bcattrib=0.5;
            else
                
                bcattrib=0.5;
            end      

        %------------------------------------------------------------------
        %Example 7.2:  This example have two material. The first is 
        %isotropic and the second one orthotropic. Dirichlet's boundary 
        %condition obtained from analitical solution (Drowniou and Le 
        %Potier, 2011). Example 4.2.2 (Eq. 53)
         
        case 7.2
            %Attribute the simple boundary condition
            bcattrib = coordmid(1);

        %------------------------------------------------------------------
        %Solve cases with SOURCE TERM
        %------------------------------------------------------------------
        %Example 11: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned by 
        %HYMAN et al, 1997)
        case 11
            %Calculate the boundary value
            bcattrib = exp(coordmid(1)*coordmid(2));
            
        %------------------------------------------------------------------
        %Example 12.1: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 1)
        case 12.1
            %Definition of parameters
            alpha = 1;
            %Calculate the pressure field
            bcattrib = calcCase12(coordmid,alpha);
                
        %------------------------------------------------------------------
        %Example 12.2: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 10)
        case 12.2
            %Definition of parameters
            alpha = 10;
            %Calculate the pressure field
            bcattrib = calcCase12(coordmid,alpha);

        %------------------------------------------------------------------
        %Example 12.3: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 100)
        case 12.3
            %Definition of parameters
            alpha = 100;
            %Calculate the pressure field
            bcattrib = calcCase12(coordmid,alpha);

        %------------------------------------------------------------------
        %Example 12.4: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 1)
        case 12.4
            %Definition of parameters
            alpha = 1000;
            %Calculate the pressure field
            bcattrib = calcCase12(coordmid,alpha);

        %------------------------------------------------------------------
        %Example 13: proble with domain orthotropic (10:1) and homogen, 
        %Dirichlet boundary condition. Analitical solution obtained from 
        %Chen et al., 2008 (first example, pp 1712, non-structured mesh) 
        case 13
            bcattrib = (cos(pi*coordmid(1))*cos(pi*coordmid(2))) + 2;

        %------------------------------------------------------------------
        %Example 14.1: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.1 theirs (mild anisotropy), pp. 2
        case 14.1
            bcattrib = 16*coordmid(1)*(1 - coordmid(1))*coordmid(2)*(1 - ...
                coordmid(2));

        %------------------------------------------------------------------
        %Example 14.2: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.2 theirs (mild anisotropy), pp. 3
        case 14.2
            bcattrib = sin((1 - coordmid(1))*(1 - coordmid(2))) + ...
                ((1 - coordmid(1))^3)*((1 - coordmid(2))^2);

        %------------------------------------------------------------------
        %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %The parameter "delta" is 10
        case 15.1
            bcattrib = sin(pi*coordmid(1))*sin(pi*coordmid(2));

        %------------------------------------------------------------------
        %Example 15.2: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %The parameter "delta" is 100
        case 15.2
            bcattrib = sin(pi*coordmid(1))*sin(pi*coordmid(2));

        %------------------------------------------------------------------
        %Example 15.3: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
        case 15.3
            bcattrib = sin(pi*coordmid(1))*sin(pi*coordmid(2));
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % adicionado
        % Este problema foi adaptado de Sheng e Yuan 2016
        case 15.5
            bcattrib = sin(pi*coordmid(1))*sin(pi*coordmid(2));
		case 15.6
		% Este problema foi adaptado de Gao and Wu 2010		
		bcattrib=0.5*((sin((1-coordmid(1))*(1-coordmid(2)))/(sin(1)))+(1-coordmid(1))^3*(1-coordmid(2))^2);
        % problema foi adpatado do artigo Gao e Wu 2015 e o artigo de Terekhov
    % 2016
%         case 15.7
%             
%         bcattrib=0; 
       % Probelam FRIIS EDWARDS A family of MPFA finite-volume schemes with full pressure support
            % for the general tensor pressure equation on cell-centered triangular grids
        case 15.8
           bcattrib=0; 
        case 15.9
           bcattrib=16*coordmid(1)*(1-coordmid(1))*coordmid(2)*(1-coordmid(2));

        %------------------------------------------------------------------
        %Example 16:  In this example there are two material with the first 
        %isotropic and the second one orthotropic. Dirichlet's boundary 
        %condition obtained from analitical solution (Drowniou and Le 
        %Potier, 2011). Example 4.2.1 (Eq. 51 and 52)
        case 16
            %When the domain is minor than 1/2 in the direction x
            if coordmid(1) <= 0.5
                bcattrib = cos(pi*coordmid(1))*sin(pi*coordmid(2));
            %When the domain is major than 1/2 in the direction x
            elseif coordmid(1) > 0.5
                bcattrib = 0.01*cos(pi*coordmid(1))*sin(pi*coordmid(2));
            end  %End of IF

        %------------------------------------------------------------------
        %Example 17:  Middle anisotropy. Section 3.1, Case 1. Obtained from
        %Gao and Wu (2010)
        case 17
            bcattrib = ...
                0.5*((sin((1 - coordmid(1))*(1 - coordmid(2)))/(sin(1))) + ...
                ((1 - coordmid(1))^3)*((1 - coordmid(2))^2)); 
            
        %When the boundary condition cames from "bcflag"
        otherwise
            %Attribute the boundary condition's value from "bcflag"
            bcattrib = bcflag(flagptr,2);
    end  %End of SWITCH

%--------------------------------------------------------------------------
%FUNCTION
%--------------------------------------------------------------------------
%Get several function of each case

%--------------------------------------------------------------------------
%Case 4:

function [bcattrib] = ...
                calcCase4(coordmid,numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4)
%Convert the coordinate x and y to polar coordinate (r,teta)
%Obtain the equivalence to "r"
r = norm(coordmid);
%Obtain the equivalence to "teta"
%For two first quadrants 
if coordmid(2) >= 0
    %Calculate the angle for two firt quadrants
    teta = acos(coordmid(1)/r);
%For two last quadrants 
elseif coordmid(2) < 0
    %Calculate the angle for two firt quadrants
    teta = 2*pi - acos(coordmid(1)/r);
end  %End of IF

%Chose the solution as a function of "numcase"
%"angle" equal to 2*pi/3 (Case 9.1 and 9.2)
if numcase == 1
    %Define the "angle"
    angle1 = 2*pi/3;
    angle2 = 5*pi/3;
%"angle" equal to pi/2 (Case 9.3)
elseif numcase == 2
    %Define the "angle"
    angle1 = pi/2;
    angle2 = 3*pi/2;
%"angle" equal to pi/3 (Case 9.4 and 9.5)
elseif numcase == 3
    %Define the "angle"
    angle1 = pi/3;
    angle2 = 4*pi/3;
end  %End of IF (angle)

if teta >= 0 && teta < angle1
    bcattrib = (r^alfa)*(a1*cos(alfa*teta) + ...
        b1*sin(alfa*teta));
elseif teta >= angle1 && teta < pi
    bcattrib = (r^alfa)*(a2*cos(alfa*teta) + ...
        b2*sin(alfa*teta));
elseif teta >= pi && teta < angle2
    bcattrib = (r^alfa)*(a3*cos(alfa*teta) + ...
        b3*sin(alfa*teta));
elseif teta >= angle2 && teta < 2*pi
    bcattrib = (r^alfa)*(a4*cos(alfa*teta) + ...
            b4*sin(alfa*teta));
end  %End of IF

%--------------------------------------------------------------------------
%Case 5: ???????

function [bcattrib] = calcCase5(coordmid,kapla)
%Definition of important parameters:
a = (6/pi)*atan(1/sqrt(1 + (2*kapla)));
d = cos(a*pi/3)/sin(a*pi/6);
%Convert the coordinate x and y to polar coordinate (r,teta)
%Obtain the equivalence to "r"
r = norm(coordmid);
%Obtain the equivalence to "teta"
%For two first quadrants 
if coordmid(2) >= 0
    %Calculate the angle for two firt quadrants
    teta = acos(coordmid(1)/r);
%For two last quadrants 
elseif coordmid(2) < 0
    %Calculate the angle for two firt quadrants
    teta = 2*pi - acos(coordmid(1)/r);
end  %End of IF
    
%Obtain the solution to angles minor than 2pi/3
if teta > 0 && teta <= 2*pi/3
    bcattrib = (r^a)*cos(a*(teta - (pi/3)));
elseif teta > 2*pi/3 && teta <= pi
    bcattrib = (r^a)*d*sin(a*((5*pi/6) - teta));
elseif teta > pi && teta <= 5*pi/3
    bcattrib = -(r^a)*cos(a*((teta - pi) - (pi/3)));
elseif teta > 5*pi/3 && teta <= 2*pi
    bcattrib = -(r^a)*d*sin(a*((5*pi/6) - (teta - pi)));
end  %End of IF

%--------------------------------------------------------------------------
%Case 6:

function [bcattrib] = calcCase6(coordmid,kapla)
%Definition of important parameters:
a = (3/pi)*atan(sqrt(1 + (2/kapla)));
d = cos(a*pi/3)/cos(2*a*pi/3);
%Convert the coordinate x and y to polar coordinate (r,teta)
%Obtain the equivalence to "r"
r = norm(coordmid);
%Obtain the equivalence to "teta"
%Calculate the tangent of teta
tanteta = ...
    abs(coordmid(2)/coordmid(1));
    %Calculate the angle "teta"
    teta = atan(tanteta);
    %Tratment to second quadrant 
    if coordmid(1) < 0 && coordmid(2) > 0
        teta = pi - teta;
    %Third quadrant
    elseif coordmid(1) < 0 && coordmid(2) < 0
        teta = teta + pi;
    elseif coordmid(1) > 0 && coordmid(2) < 0
        teta = 2*pi - teta;
    end  %End of IF
    
    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 && teta < 2*pi/3
        bcattrib = (r^a)*cos(a*(teta - (pi/3)));
    else
        bcattrib = (r^a)*d*cos(a*((4*pi/3) - teta));
    end  %End of IF

%--------------------------------------------------------------------------
%Case 12:

function [bcattrib] = calcCase12(coordmid,alpha)
%When the coordinate x is lower than 0
if coordmid(1) < 0
    bcattrib = (((2*sin(coordmid(2))) + cos(coordmid(2)))*...
        (alpha*coordmid(1))) + (sin(coordmid(2)));
%When the coordinate x is major than 0
elseif coordmid(1) >= 0
    bcattrib = (exp(coordmid(1)))*sin(coordmid(2));
end  %End of IF
