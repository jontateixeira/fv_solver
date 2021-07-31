function [transmvecleft,transmvecright,knownvecleft,knownvecright,...
          storeinv,Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
          bodyterm] = transmTPS(kmap,knownboundlength)
%Define global parameters:
global coord nsurn1 phflw

if phflw==1 | phflw==0
    phasekey=1;
elseif phflw==2
    phasekey=2;
end

storeinv = 0;
Bleft = 0;
Bright = 0;
bodyterm = 0;
pointedge = 0;
Fg = 0;

%Initialize "transmvecleft". This parameter stores the little matrix calc.
%by the procedure [A] + [B][D]-1[C]. This will be used to calculate the 
%velocity field. The matrix "[T] = [A] + [B][D]-1[C]" will be stored as 
%follows: each row of "T" will be stored into vector "transmvecleft", row 
%to row. All matrix will alocate "nsurn*esurn" positions in the vector.
%Obs.: ONLY THE CONTRIBUTION ON THE LEFT IS STORED (only the velocity on 
%the left is calculated).
[amountesurn,amountnsurn] = catchlength;
%Define the "transmvecleft" length
lengthstoreinv = (amountesurn')*amountnsurn;
%On the left
transmvecleft = zeros(lengthstoreinv,1);
%On the right
transmvecright = transmvecleft;
%Initialize "maptransm". It stores the pointer to "transmvec"
maptransm = zeros(size(coord,1) + 1,1);
mapinv = maptransm;
mapknownvec = maptransm;

%Initialize "knownvecleft". This vector stores the known values, initialy
%attributed to leftvector.
knownvecleft = zeros(length(nsurn1),1);
knownvecright = knownvecleft;

%Initialize "incrementrowtransm" and "incrementvec". This parameter 
%increment according the row's number the vectors "transmvecleft" and 
%"other", respectively.
incrementrowtransm = 0;
incrementrowinv = 0;
incrementvec = 0;

%--------------------------------------------------------------------------
%Swept all nodes (from 1 until nnode)
    
for inode = 1:size(coord,1)
    %Swept the amount of elements surrounding each node evaluated 
    %using the linked list "esurn2" pointing to "esurn1".
    %Each interaction region is centered in the node evaluated.

    %Obtain "esurn" and "nsurn" to each node evaluated. Use the
    %function "getsurnode"
    [esurn,nsurn] = getsurnode(inode);
        
    %----------------------------------------------------------------------
    %Call the function "calclittlematrix". This function calculate the
    %local matrix which later will be assembled into global matrix "M"

    %It gets the little matrix
    [subunknAleft,subunknAright,subauxBleft,subauxBright,subunknC,...
      subauxD,getinfoedge,Fg] = calclittlematrix(inode,kmap,phasekey,...
      esurn,nsurn,Fg);
    
    %Some parameters must be initialized accordin the matrix calculated
        
    %Initialize the little vectors ("subvector") whose will fill the 
    %"mvector" in each node evaluated. ("subvector" just will have 
    %nonzero values if either Dirichlet or Newmann boundary condition 
    %are applied over half edges). Is possible have "subvector" till in
    %situations which "amountunknow" is null. 
    subvectorleft = zeros(length(nsurn),1);
    %The same occur to "subvectoright"
    subvectoright = subvectorleft;
    %"knowntherm" receives the therm of "subauxBleft" multiplied by
    %known pressure, This occur when a Dirichlet boundary condition is
    %appliesd to surface (or half edge in case of two-dimensional
    %domain). The order of vector "knowntherm" has the same order of 
    %"subvectorleft"
    knowntherm = zeros(length(nsurn),1);
    %"knownflow" receives all known flow rate (non null Newmann 
    %boundary condition) multiplied by inv(subauxD). See below
    knownflow = zeros(length(nsurn),1);
    %"flowvalue" receives the value of each Newmann boundary condition
    flowvalue = zeros(length(nsurn),1);
        
    %----------------------------------------------------------------------
    %There is several types of boundary conditions application:
    %Below each one of them is defined:
       
    %1. Least one half edge is under either Dirichlet or Newmann 
    %boundary condition. The procedure necessary to implementation of 
    %this kind of boundary condition will be developed
    if any(getinfoedge(:,1))  
        %Abount boundary condition type:
        %Verify which flag is of Dirichlet. This parameter takes
        %the row of "getinfoedge" with this feature.
        dirichpointer = ...
            logical(getinfoedge(:,1) > 0 & getinfoedge(:,1) < 200);
        %Verify which flag is of Neumann. This parameter takes
        %the row of "getinfoedge" with this feature.
        neumannpointer = ...
            logical(getinfoedge(:,1) > 200 & getinfoedge(:,1) < 300);
                        
        %Once there is a Dirichlet boundary condiion, a set of
        %operations will be done
        if any(dirichpointer)
            %Change the auxilary matrix "subauxD" according to
            %auxilary pressure with Dirichlet boundary condition.
                
            %Attribute zero for "subauxD" rows associated with 
            %auxlary pressure with Dirichlet boundary condition.
            subauxD(dirichpointer,:) = 0;
            %Attribute zero for "subunknC" rows associated with 
            %auxlary pressure with Dirichlet boundary condition.
            subunknC(dirichpointer,:) = 0;
            %"knowntherm" receives the Dirichlet boundary condition 
            %values in its respective position
            knowntherm(dirichpointer) = getinfoedge(dirichpointer,2);
            %After to turn null the row which corresponds to Dirichlet 
            %boundary condition, attribute 1 to diagonal associated 
            %with auxiliary pressure with Dirichlet boundary condit.
            subauxD(logical(diag(dirichpointer))) = 1;
        end  %End of internal IF (Dirihlet)
        %Once there is a Newmann boundary condiion, a set of
        %operations will be done
        if any(neumannpointer)
            %"knownflow" receives the Newmann boundary condition 
            %values in its respective position
            flowvalue(neumannpointer) = ...
                (getinfoedge(neumannpointer,2).*...
                getinfoedge(neumannpointer,3))/knownboundlength;
        end  %End of internal IF (Neumann)
                        
        %Obtain the inverse of "subauxD" matrix
        invmatrix = subauxD\eye(size(subauxD,1));
        %Once the changed matrix is ready, it is inverted, after 
        %that, it will be multiplied by "subauxC".
        %The result of that is attibuted to "subunknaux".
        subunknaux = invmatrix*subunknC;
        %"knowntherm" receives influence of "inv(subauxD)". Dirichlet
        knowntherm = invmatrix*knowntherm;
        %"knownflow" receives influence of "inv(subauxD)". Neumann
        knownflow = invmatrix*flowvalue;
            
        %Fill the matrix ([A] + [B]{P-}), where {p-} = [D]-1*[C]. 
        %Later this matrix is associated to global matrix "M" 
        %Fill "subunknleft"
        subunknleft = subunknAleft + (subauxBleft*subunknaux);
        %Fill "subunknright"
        subunknright = subunknAright + (subauxBright*subunknaux);
        %Fill the vector which is later associated to global vector, 
        %"mvector"
        %Fill "subvectorleft", contribution of left swing (Dirichlet)
        subvectorleft = subvectorleft + (subauxBleft*knowntherm);
        %Fill "subvectorleft", contribution of left swing (Newmann)
        subvectorleft = subvectorleft + (subauxBleft*knownflow);

        %Fill "subvectoright", contribution of right swing (Dirichlet)
        subvectoright = subvectoright + (subauxBright*knowntherm);
        %Fill "subvectoright", contribution of left swing (Newmann)
        subvectoright = subvectoright + (subauxBright*knownflow);
            
    %2. Iteraction Region without boundary condition.
    else
        %Obtain the inverse of "subauxD"
        invmatrix = subauxD\eye(size(subauxD,1));
        %Building "subunknleft" (subunknleft = A + B*(D^-1)*C, or 
        %TRANSMISIBILITY)
        subunknleft = subunknAleft + (subauxBleft*(invmatrix*subunknC));
        %Building "subunknright" (subunknright = A + B*(D^-1)*C, or 
        %TRANSMISIBILITY)
        subunknright = subunknAright + (subauxBright*(invmatrix*subunknC));
    end  %End of IF (boundary or internal edges)

    %----------------------------------------------------------------------
    %It stores the matrices (boundary or inside domain):
            
    %It stores the matrices:
    %Store in "transmvecleft" the aux. matrix "[A] + [B][D]-1[C]"
    transmvecleft(incrementrowtransm + 1:incrementrowtransm + ...
        length(esurn)*length(nsurn)) = subunknleft';
    %Store in "transmvecright" the aux. matrix "[A] + [B][D]-1[C]"
    transmvecright(incrementrowtransm + 1:incrementrowtransm + ...
        length(esurn)*length(nsurn)) = subunknright';

    %It stores vectors:
    %Store in "storevec" the known vector "subvectorleft"
    knownvecleft(incrementvec + 1:incrementvec + ...
        length(nsurn),1) = subvectorleft;
    %Store in "storevec" the known vector "subvectorright"
    knownvecright(incrementvec + 1:incrementvec + ...
        length(nsurn),1) = subvectoright;

    %----------------------------------------------------------------------
    %Increment "incrementrowinv", "incrementrowtransm" and 
    %"incrementvec" with row's length of little matrix stored
    incrementrowinv = incrementrowinv + length(nsurn)^2;
    incrementrowtransm = incrementrowtransm + ...
        (length(esurn)*length(nsurn));
    incrementvec = incrementvec + length(nsurn);

    %Fill "maptransm" and "mapknownvec" with "incrementrowtransm" and 
    %"incrementvec" positions
    mapinv(inode + 1) = incrementrowinv;
    maptransm(inode + 1) = incrementrowtransm;
    mapknownvec(inode + 1) = incrementvec; 
end  %End of external FOR (Swept each node - calculate the pressure)

%It clears the vectors "ibedg" and "iinedg".
clear ibedg iinedg;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "midedge"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordmidedge] = midedge(inode,nodeval)
%Define global parameters:
global coord;

%Obtain the vector "coordmidedge" using a mean between the two vectors 
%which constitute the edge evaluated 
coordmidedge = 0.5*(coord(inode,:) + coord(nodeval,:));

%--------------------------------------------------------------------------
%FUNCTION "definepoints"
%--------------------------------------------------------------------------

function [xcol,ycol,xneta,yneta,xzeta,yzeta] = definepoints(elemeval,...
    inode,netapoint,zetapoint)
%Definition of global parameters
global centelem;

%Coordinate "x" of colocation point (left element)
xcol = centelem(elemeval,1);
%Coordinate "y" of colocation point (left element)
ycol = centelem(elemeval,2);
%Coordinate of midpoint "w" (netanode)
%Calculate the midedge coordinate using the function "midedge"
netacoord = midedge(inode,netapoint);
%Obtain the coordinate "x"
xneta = netacoord(1);
%Obtain the coordinate "y"
yneta = netacoord(2);
%Coordinate of midpoint "s" (zetanode)
%Calculate the midedge coordinate using the function "midedge"
zetacoord = midedge(inode,zetapoint);
%Obtain the coordinate "x"
xzeta = zetacoord(1);
%Obtain the coordinate "y"
yzeta = zetacoord(2);

%--------------------------------------------------------------------------
%FUNCTION "calcinvj"
%--------------------------------------------------------------------------

%The function "calcinvj" calculate the Jacobian's determinant. 
function [invj] = calcinvj(inode,elemeval,zetanode,netanode)
%Call "definepoints" function
[xcol,ycol,xneta,yneta,xzeta,yzeta] = definepoints(elemeval,inode,...
    netanode,zetanode);

%Initialize the jacobian matrix
j = zeros(2);
%Build the jacobian matrix
%dx/dzeta
j(1,1) = (xzeta - xcol);
%dy/dzeta
j(1,2) = (yzeta - ycol);
%dx/dneta
j(2,1) = (xneta - xcol);
%dy/dneta
j(2,2) = (yneta - ycol);

%Calculate the inverse of jacobian matrix (left element)
invj = inv(j);

%--------------------------------------------------------------------------
%FUNCTION "calcgradP"
%--------------------------------------------------------------------------

%Calculate the gradient of pressure to both right and left elements 
%This returns a matrix [2x4] for each element evaluated
function [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
    inode,zetanode,netanode,edgekey)
%Initialize "gradPleft"
gradPleft = zeros(3);
%Initialize "gradPright"
gradPright = 0;

%Define "gradmap". This parameter is the matrix which denots grad(P) in
%natural coordinate (zeta,neta)
gradmap = zeros(2,3);
%Fill "gradmap"
gradmap(1,:) = [-1 1 0];
gradmap(2,:) = [-1 0 1];

%Define "gradPleft" multiplying the inverse of jacobian matrix per 
%"gradmap". Obtain the inverse of jacobian
invj = calcinvj(inode,leftelemeval,zetanode,netanode(1));
%Fill "gradPleft"
gradPleft(1:2,:) = invj*gradmap;

%When the edge evaluated is inside domain, there is necessity of evaluate
%one element to left and other element to right. In this case:
if edgekey == 1
    %Initialize "gradPright"
    gradPright = zeros(3);
    %Define "gradPright" multiplying the inverse of jacobian matrix per 
    %"gradmap".
    %Obtain the inverse of jacobian (right element)
    invj = calcinvj(inode,rightelemeval,zetanode,netanode(2));
    %Fill "gradPright"
    gradPright(1:2,:) = invj*gradmap;
end  %End of IF

%--------------------------------------------------------------------------
%Function "calclittlematrix"
%--------------------------------------------------------------------------

function [subunknAleft,subunknAright,subauxBleft,subauxBright,subunknC,...
          subauxD,getinfoedge,Fg] = calclittlematrix(inode,kmap,phasekey,...
          esurn,nsurn,Fg)
%Define global parameters:
global elem bedge inedge normals bcflag g rowposit nsurn2; 

%Initializa "getinfoedge". This parameter gives informations about
%boundary condition of each half edge which surrounding the node evaluated.
%If the edge isn't over the boundary "getinfoedge" will receive zero in its 
%first column. Otherwise, this parameter receives the code associate with 
%each boundary condition stabeleted over the edge. In the second column is 
%usualy stored the value of boundary condition, obtained by "bcflag" table. 
%In case of there ain't boundary condition (edge evaluated inside domain)
%the second column's value keeping zero too. Finaly the third column stores
%the normal length (case need).
getinfoedge = zeros(length(nsurn),3);

%The matrix "subunknAleft" and "subauxBleft" are used in the equality 
%f = AP^ + BP-. The flow rate in this equation is obtained by 
%contribution of elements to left of the half edge evaluated.
        
%Initialize the matrix "subunknAleft" (matrix A) and "subauxBleft" 
%(matrix B) which countain the coefitients associated respectivaly to 
%unknown node and auxilary nodes.
%The order of matrix "subunknAleft" is referred to among of nodes which are 
%not under Dirichlet boundary condition x amount of elements which concorre 
%to node evaluated
subunknAleft = zeros(length(nsurn),length(esurn));
%"subunknAright" has the same order of "subunknAleft" and account the 
%contribution of flow rate over the elements to right of the half edge 
%evaluated.
subunknAright = subunknAleft;
%"subauxBleft" has the order equal to number of edges where the flow rate 
%must be calculated.
subauxBleft = zeros(length(nsurn));
%"subauxBright" is the same
subauxBright = subauxBleft;
%The matrix "subunknC" and "subauxD" are used in the equality 
%f = CP^ = DP-. The matrix "subunknC" and "subauxD" are initialized.
%Considering which the referred matrix has contribution of both 
%left and right elements is not necessary have more than one matrix
%to this finality.
subunknC = zeros(length(nsurn),length(esurn));
subauxD = zeros(length(nsurn));

%Initialize bodyterms. It stores gravity contribution for local algebric 
%system 
bodytermleft = zeros(length(nsurn),1);
bodytermright = bodytermleft;
edgemap = bodytermleft;

%Initialize "Krock"
kleft = zeros(3);
kright = kleft;

%Get the position of row for the interaction region eval.
rowposbyir = rowposit(nsurn2(inode) + 1:nsurn2(inode + 1));

%Swept AGAIN all half edges (or surfaces) surround each node evaluated
for insurn = 1:length(nsurn)
    %Node which surround the node evaluated. This is 
    %obtained from nsurn1 with pointer from nsurn2
    inodesurn = nsurn(insurn);
    %Define the "bedge" size
    bedgesize = size(bedge,1);

    %The HALF-EDGE BELONGS to BOUNDARY edges 
    if rowposbyir(insurn) <= bedgesize
        %Get "bedgerow"
        bedgerow = rowposbyir(insurn);
        %In the edges over boundary there is only "leftelemeval" 
        leftelemeval = bedge(bedgerow,3);
        
        %"zetanode" receives the auxilary node which is not 
        %"inodesurn". Case the edge is inner domain, a vector [1x2] 
        %gives this parameter. the first is the left "zetanode" and the 
        %second is the right "zetanode"  
        netanode = calcnetanode(inodesurn,nsurn,leftelemeval,0,0);
        %"netanode" receives "inodesurn" due this node define the edge
        %evaluated.
        zetanode = inodesurn;
            
        %"calcnflow" is a function which calculate the normal to 
        %half surfaces (or half edge) of the elements to left 
        nflowleft = 0.5*normals(bedgerow,:);

        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [2x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,] = calcgradP(leftelemeval,0,inode,zetanode,...
           netanode,0);        
                
        %"getinfoedge" gives information about the boundary
        %condition. Its first column is referred to flag imposed to
        %edge. The second one is referred to Boundary condition's 
        %value attributed to it.
        getinfoedge(insurn,1) = bedge(bedgerow,5);

        %Store in the second column of "getinfoedge" the value of
        %boundary condition. 
        %Define "flagpointer"
        flagpointer = logical(bcflag(:,1) == bedge(bedgerow,5));
        %It Catches the boundary condition value
        getinfoedge(insurn,2) = PLUG_bcfunction([inode inodesurn],...
            flagpointer);
          
        %Attribute to third "getinfoedge"'s column the "normals"'s length 
        getinfoedge(insurn,3) = norm(nflowleft);

        %------------------------------------------------------------------
        %Fill matrix of permeability as a function of element 
        %to left of half edge evaluated.

        %"kfeaturelef" follows the same procedure which internal nodes. 
        kfeaturelef = elem(leftelemeval,5);
                
        %Obtain the permeability of element on the left of half edge.
        kleft(1:2,1:2) = ...
            [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];

        %------------------------------------------------------------------
        %BODY FORCE PARAMETERS:

        %In two-phase case ("phasekey" == 2)
        if norm(g) ~= 0
            %Left contribution
            prodkg = kleft*g;
            %Attibute "dotkgn" to "bodytermleft"
            bodytermleft(insurn) = dot(prodkg',nflowleft);
            %Attribute to "edgemap" the edge number
            edgemap(insurn) = bedgerow;
            
            %Body force parameter used in hyperbolic equation (Fg)
            if phasekey == 2
                %Attribute the dot product K*g.n to "Fg"
                Fg(leftelemeval,:) = Fg(leftelemeval,:) + ...
                    prodkg(1:2)'*(any(Fg(leftelemeval,:)) == 0);
            end  %End of IF (hyperbolic equation contribution)
        end  %End of IF (body force in hyperbolic equation)

        %------------------------------------------------------------------
        %Fill both the matrix "subunknAleft" and the matrix "subauxBleft"
            
        %Calculate "k.nablaP" to left element
        knablap = kleft*gradPleft;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAleft" and is associated to COLOCATION point
        colocterm = -nflowleft*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %The matrix "subunknA" is filled: 
        subunknAleft(insurn,findcol) = subunknAleft(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to ZETA point
        zetaterm = -nflowleft*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to NETA point
        netaterm = -nflowleft*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode);
        %Second auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            netaterm;

        %------------------------------------------------------------------
        %Fill both the matrix "subunknC" and the matrix "subauxD" 

        %Contribution of element to LEFT of half edge evaluated: 
        %COLOCATION matrix
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAleft(insurn,:);
        %Contribution of element to LEFT of half edge evaluated: 
        %AUXILARY matrix
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBleft(insurn,:);
        
    %The HALF-EDGE BELONGS to INERNAL edges        
    else
        %Get "bedgerow"
        inedgerow = rowposbyir(insurn) - bedgesize;
        
        %"leftelemeval" is the element evaluated to left of 
        %each half edge evaluated
        leftelemeval = inedge(inedgerow,3);
        %"rightelemeval" is the element evaluated to right of 
        %each half edge evaluated
        rightelemeval = inedge(inedgerow,4);
                
        %"zetanode" receives the auxilary node which is not 
        %"inodesurn". Case the edge is inner domain, a vector [1x2] 
        %gives this parameter. the first is the left "zetanode" and the 
        %second is the right "zetanode"  
        netanode = calcnetanode(inodesurn,nsurn,leftelemeval,...
            rightelemeval,1);
        %"netanode" receives "inodesurn" due this node define the edge
        %evaluated.
        zetanode = inodesurn;
            
        %Get the normal vector to half-edge evaluated (left element).
        nflowleft = 0.5*normals(bedgesize + inedgerow,:);
        %Get the normal vector to half-edge evaluated (right element).
        nflowright = -nflowleft;
            
        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [2x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
            inode,zetanode,netanode,1);
            
        %------------------------------------------------------------------
        %Obtain the tensor permeability as a function of element 
        %to left and right of half edge evaluated

        %"kfeaturelef" associates the values of permeability 
        %tensor with the element to left of half edge evaluated. 
        %Depending of domain, each element may have a 
        %permeability map. 
        kfeaturelef = elem(leftelemeval,5);
        %"kfeaturerig" associates the value of permeability 
        %tensor with the element to right of half edge 
        %evaluated. Depending of domain, each element may have  
        %a permeability map. 
        kfeaturerig = elem(rightelemeval,5);
            
        %Obtain the permeability of element on the left of half edge.
        kleft(1:2,1:2) = ...
            [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];
            
        %Obtain the permeability of element on the left of half edge.
        kright(1:2,1:2) = ...
            [kmap(kfeaturerig,2:3); kmap(kfeaturerig,4:5)];

        %------------------------------------------------------------------
        %BODY FORCE PARAMETERS:

        %In two-phase case ("phasekey" == 2)
        if norm(g) ~= 0
            %Left contribution
            prodkgleft = kleft*g;
            %Attibute "dotkgn" to "bodytermleft"
            bodytermleft(insurn) = dot(prodkgleft',nflowleft); 
            %Right contribution
            prodkgright = kright*g;
            %Attibute "dotkgn" to "bodytermright"
            bodytermright(insurn) = dot(prodkgright',nflowright); 
            %Attribute to "edgemap" the edge number
            edgemap(insurn) = bedgesize + inedgerow;

            %Body force parameter used in hyperbolic equation (Fg)
            if phasekey == 2
                %Attribute the dot product K*g.n to "Fg"
                Fg(leftelemeval,:) = Fg(leftelemeval,:) + ...
                    prodkgleft(1:2)'*(any(Fg(leftelemeval,:)) == 0);
                %Attribute the dot product K*g.n to "Fg"
                Fg(rightelemeval,:) = Fg(rightelemeval,:) + ...
                    prodkgright(1:2)'*(any(Fg(rightelemeval,:)) == 0);
            end  %End of IF
        end  %End of IF (body force in hyperbolic equation)

        %------------------------------------------------------------------
        %Fill both the matrix "subunknAleft" and the matrix 
        %"subauxBleft"
        
        %Calculate only the flows in the LEFT of each half edge 
        %evaluated in order honer the equality f = AP^ + BP-. The
        %procedure follows.
            
        %Calculate "k.nablaP" to left element
        knablap = kleft*gradPleft;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAleft" and is associated to COLOCATION point
        colocterm = -nflowleft*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %The matrix "subunknA" is filled: 
        subunknAleft(insurn,findcol) = subunknAleft(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to ZETA point
        zetaterm = -nflowleft*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to NETA point
        netaterm = -nflowleft*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(1));
        %Second auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            netaterm;
            
        %------------------------------------------------------------------
        %Fill both the matrix "subunknAright" and the matrix 
        %"subauxBright"

        %Calculate only the flows in the RIGHT of each half edge 
        %evaluated in order honer the equality f = AP^ + BP-. 
            
        %Calculate "k.nablaP" to RIGHT element
        knablap = kright*gradPright;

        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAright" and is associated to COLOCATION point
        colocterm = -nflowright*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"rightelemeval". "findcol" receives the index.
        findcol = logical(esurn == rightelemeval);
        %The matrix "subunknA" is filled: 
        subunknAright(insurn,findcol) = subunknAright(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBright" and is associated to ZETA point
        zetaterm = -nflowright*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"rightelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge zeta)
        subauxBright(insurn,findcol) = subauxBright(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBright" and is associated to NETA point
        netaterm = -nflowright*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"rightelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(2));
        %Second auxilary node (half edge neta)
        subauxBright(insurn,findcol) = subauxBright(insurn,findcol) + ...
            netaterm;

        %------------------------------------------------------------------
        %Fill both the matrix "subunknC" and the matrix "subauxD"

        %Contribution of element to LEFT of half edge evaluated: 
        %COLOCATION matrix. The value attributed to "subunknC" is
        %the same attibuted to "subunknAleft" 
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAleft(insurn,:);

        %Contribution of element to RIGHT of half edge evaluated:
        %COLOCATION matrix. The value attributed to "subunknC" is
        %the same attibuted to "subunknAright"
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAright(insurn,:);

        %Contribution of element to LEFT of half edge evaluated: 
        %AUXILARY matrix. This value is the same (with sign changed) 
        %of that attributed to "subauxBleft" 
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBleft(insurn,:);
 
        %Contribution of element to RIGHT of half edge evaluated: 
        %AUXILARY matrix. This value is the same (with sign changed) 
        %of that attributed to "subauxBleft" 
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBright(insurn,:);
            
    end  %End of IF (to be or not to be over the boundary) 
end  %End of first FOR (swept the nodes surrounding the node eval.)
    
%--------------------------------------------------------------------------
%Function "catchlength"
%--------------------------------------------------------------------------

function [amountesurn,amountnsurn] = catchlength
%Define global parameters:
global esurn2 nsurn2;

%Initialize "amountesurn" and "amountnsurn"
amountesurn = zeros(length(esurn2) - 1,1);
amountnsurn = zeros(length(nsurn2) - 1,1);
%It catches the amount of elements surrounding each node 
i = 1:length(esurn2) - 1;
amountesurn(i) = esurn2(i + 1) - esurn2(i);
%It catches the amount of half-edges surrounding each node 
j = 1:length(nsurn2) - 1;
amountnsurn(j) = nsurn2(j + 1) - nsurn2(j);

