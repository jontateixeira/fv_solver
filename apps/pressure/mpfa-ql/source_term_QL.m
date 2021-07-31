function [ source ] = source_term_QL( centelem, source )
%
global elem elemarea

for i=1:size(elem,1)
    x=centelem(i,1); y=centelem(i,1);
    source(i,1) = (-4.9342*(x^2-y^2)*cos(pi*y)*cos(pi*x)+...
                  (pi^2)*(2+2.9998*x^2+2.9998*y^2)*sin(pi*y)*sin(pi*x)-...
                  (1.57061*x+6.70353*y)*cos(pi*y)*sin(pi*x)+...
                  (-1.57061*y+6.70353*x)*cos(pi*x)*sin(pi*y))*(elemarea(i));
end

end

