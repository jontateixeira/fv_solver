function [knownboundlength] = getknownboundlength(klb)
%Define global parameters:
global bcflag normals;

%Initialize "knownboundlength"
knownboundlength = 0;

%It points to a Neumann Boundary Condition 
pointnonnullflag = (bcflag(:,1) > 200 & bcflag(:,2) > 0);
flagref = bcflag(logical(pointnonnullflag),1);
%Choose according "pointnonnullflag"
%There is a Neumann Boundary Condition
if any(klb)
    %Swept the amount of different flags
    for iflag = 1:length(flagref)
        %Swept all edges and get its lengths
        for i = 1:length(klb)
            iedge = klb(i);
            %Attribute the length of each edge
            knownboundlength = knownboundlength + norm(normals(iedge,1:2));
        end  %End of FOR
    end  %End of FOR
%There is NO a Neumann Boundary Condition
else
    knownboundlength = 1;
end  %End of IF
    

