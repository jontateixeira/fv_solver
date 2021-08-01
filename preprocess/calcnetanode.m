function [netanode] = calcnetanode(inodesurn,nsurn,leftelemeval,rightelemeval,edgekey)
%Define global parameters
global elem;

%"interselemnode" is a parameter which account the intersection among the 
%nodes which constitute the element and the node which surround the node 
%evaluated. 

%Initialize "nsurnptr" and "vertelem"
vertelem = elem(leftelemeval,1:4);
nsurnptr = zeros(length(nsurn),1);
%Swept "nsurn"
for i = 1:length(nsurn)
    pointinters = logical(nsurn(i) == vertelem);
    nsurnptr(i) = sum(pointinters);
end  %End of FOR
interselemnode = nsurn(logical(nsurnptr));

%Obtain a logical reference of "zetanode"
netanodepointer = logical(interselemnode ~= inodesurn);
%Obtain the "zetanode" value
netanode = interselemnode(netanodepointer);

%Case there is a RIGHT element
if edgekey == 1
    %Obtain "interselemnode" to element to right
    %Initialize "nsurnptr" and "vertelem"
    vertelem = elem(rightelemeval,1:4);
    nsurnptr = zeros(length(nsurn),1);
    %Swept "nsurn"
    for i = 1:length(nsurn)
        pointinters = logical(nsurn(i) == vertelem);
        nsurnptr(i) = sum(pointinters);
    end  %End of FOR
    interselemnode = nsurn(logical(nsurnptr));
    
    %Obtain a logical reference of "zetanode"
    netanodepointer = logical(interselemnode ~= inodesurn);
    %Obtain the "zetanode" value
    netanode(2) = interselemnode(netanodepointer);
end  %End of IF
