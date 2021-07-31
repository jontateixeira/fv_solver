function [mtxcoeff,vecoeff,esurn] = gettransmcoeff(transmvecleft,...
    knownvecleft,storeinv,B,flagside,mapinv,maptransm,mapknownvec,...
    pointedge,bodyterm,vrtx1,vrtx2,phasekey,row,posvertex,mobility,...
    bedgesize,inedgesize)
%Get the mobility through edge
totalmobility = getotalmobility(mobility,phasekey,row,posvertex,bedgesize,...
    inedgesize);
%Get the terms of matrix [A] + [B][G]-1[H] (left)
mtxleft = transmvecleft(maptransm(vrtx1) + 1:maptransm(vrtx1 + 1));
    
%Get terms of known vector:
%Left
vecleft = knownvecleft(mapknownvec(vrtx1) + 1:mapknownvec(vrtx1 + 1));
        
%Get the elements and the nodes surrounding the vertex 1
[esurn,nsurn] = getsurnode(vrtx1);
%Find the "vertices(2)" position in "nsurn"
nsurnpos = find(nsurn == vrtx2);
%Get the nsurn'th row in "mtxleft" and "vecleft"
%Define the position.
finalpos = nsurnpos*length(esurn);
initpos = finalpos - length(esurn) + 1;
%Rebuild matrices:
mtxcoeff = totalmobility*mtxleft(initpos:finalpos);
%Rebuild vectors:
vecoeff = totalmobility*vecleft(nsurnpos);

%When there exists GRAVITY
if any(storeinv)
    %Initialize "auxcoeff"
    auxcoeff = 0;
    
    %It catches the number of rows associated to vertex evaluated.
    getnsurn = pointedge(mapknownvec(vrtx1) + 1:mapknownvec(vrtx1 + 1));
    %It catches the body force contribution on the left and on the right
    localbodyterm = ...
        bodyterm(mapknownvec(vrtx1) + 1:mapknownvec(vrtx1 + 1),:);
    
    %Calculate the local mobility vector. It has the "getnsurn" length.
    localmobvec = getotalmobility(mobility,phasekey,row,posvertex,...
        bedgesize,inedgesize);
    %Calculate a average density. It has the same "getnsurn" length.
    averdens = getdensity(mobility,localmobvec,phasekey,getnsurn,...
        posvertex,bedgesize,inedgesize);
    
    %Update "localbodyterm" with two-phase contribution (left hand side).
    localbodyterm(:,1) = averdens.*(localmobvec.*localbodyterm(:,1));
    %Update "localbodyterm" with two-phase contribution (right hand side).
    localbodyterm(:,2) = averdens.*(localmobvec.*localbodyterm(:,2));
    
    %Recovery the inv(D) matrix.
    recinvd = storeinv(mapinv(vrtx1) + 1:mapinv(vrtx1 + 1));
    %Put "recinvd" in matrix form
    invd = reshape(recinvd,length(getnsurn),length(getnsurn))';
    
    %Put the mobility contribution in "invd".
    for j = 1:length(localmobvec)
        invd(:,j) = invd(:,j)./localmobvec(j);
    end  %End of FOR
    
    %Get the product invd*sum(localbodyterm)
    G = invd*(sum(localbodyterm,2));
    
    %At last we take the corrected row of "B" matrix
    %Recovery "B"
    recoveryb = B(mapinv(vrtx1) + 1:mapinv(vrtx1 + 1));
    
    %Define the position.
    finalpos = nsurnpos*length(nsurn);
    initpos = finalpos - length(nsurn) + 1;
    %mob*[Brow]*[invD]*g
    bodycoeff = (totalmobility*recoveryb(initpos:finalpos)')*G;

    %Evaluate the side contribution (for left hand side)
    auxcoeff = auxcoeff + localbodyterm(nsurnpos,1)*(flagside == 1);
    %Evaluate the side contribution (for right hand side)
    auxcoeff = auxcoeff + localbodyterm(nsurnpos,2)*(flagside == 2);
  
    %Update "vecoeff"
    vecoeff = vecoeff + (bodycoeff + auxcoeff); 
end  %End of IF (get the average density)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getmobility"
%--------------------------------------------------------------------------

function [totalmobility] = getotalmobility(mobility,phasekey,row,posvertex,...
    bedgesize,inedgesize)
%Define global parameters
global visc;
                
%Initialize "totalmobility"
totalmobility = zeros(length(row),1);

%Define "mobility"
%In one-phase case ("phasekey" == 1)
if phasekey == 1
    totalmobility(:) = 1/visc(1);
%In two-phase case ("phasekey" == 2) with NO Kozdon Multidimensional Scheme
elseif phasekey == 2 && size(mobility,1) == bedgesize + inedgesize
    %Define "totalmobility"
    totalmobility = sum(mobility(row,:));  
%In two-phase case ("phasekey" == 2) with Kozdon MultiD Scheme
elseif phasekey == 2 && size(mobility,1) > bedgesize + inedgesize
    %Define aboolean parameter:
    booleanvtx = (posvertex == 1);
    %Water mobility when "inode" is lower than "inodesurn" 
    %("booleanvtx" = 1) or when "inode" is bigger than "inodesurn" 
    %("booleanvtx" = 0)
    watermobility = mobility(2*row - 1,1)*booleanvtx + ...
        mobility(2*row,1)*(1 - booleanvtx);

    %Oil mobility when "inode" is lower than "inodesurn" (1) or
    %when "inode" is bigger than "inodesurn" (0)
    oilmobility = mobility(2*row - 1,2)*booleanvtx + ...
        mobility(2*row,2)*(1 - booleanvtx);

    %Define "totalambda" for half-edge evaluated
    totalmobility = watermobility + oilmobility;
end  %End of IF

%--------------------------------------------------------------------------
%Function "getdensity"
%--------------------------------------------------------------------------

function [averdens] = getdensity(mobility,localmobvec,phasekey,row,...
    posvertex,bedgesize,inedgesize)
%Define global parameters:
global dens;

%Initialize "averdens"
averdens = zeros(length(row),1);

%Define the average density:
%In one-phase case ("phasekey" == 1)
if phasekey == 1
    averdens(:) = dens(1);
%In two-phase case ("phasekey" == 2) with NO Kozdon Multidimensional Scheme
elseif phasekey == 2 && size(mobility,1) == bedgesize + inedgesize
    %Define "mobility"
    averdens = (mobility(row,:)*dens')./localmobvec;
%In two-phase case ("phasekey" == 2) with Kozdon MultiD Scheme
%???????????????????????????????????????????????????????????????????????
elseif phasekey == 2 && size(mobility,1) > bedgesize + inedgesize
    %Initialize phase mobilities
    watermobility = 0;
    oilmobility = 0;
    %Water mobility (when "inode" is lower than "inodesurn")
    watermobility = watermobility + ...
        mobility(2*row - 1,1)*(posvertex == 1);
    %Water mobility (when "inode" is bigger than "inodesurn")
    watermobility = watermobility + ...
        mobility(2*row,1)*(posvertex == 2);

    %Oil mobility (when "inode" is lower than "inodesurn")
    oilmobility = oilmobility + ...
        mobility(2*row - 1,2)*(posvertex == 1);
    %Oil mobility (when "inode" is bigger than "inodesurn")
    oilmobility = oilmobility + ...
        mobility(2*row,2)*(posvertex == 2);

    %Define "totalambda" for half-edge evaluated
    averdens = watermobility*dens(1) + oilmobility*dens(2);
end  %End of IF
