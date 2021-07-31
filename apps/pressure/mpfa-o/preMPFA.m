function [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
    Bleft,Bright,mapinv,maptransm,mapknownvec,pointedge,bodyterm] = preMPFA(kmap,klb)

%Obtain the coordinate of both CENTER and AUXILARY

%Get the length of the edge with non-null Neumann Boundary Condition.
knownboundlength = getknownboundlength(klb);

%--------------------------------------------------------------------------
[transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmTPS(kmap,1,knownboundlength);

