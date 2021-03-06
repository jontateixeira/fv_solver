%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%--------------------------------------------------------------------------

function [coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
        normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
        satlimit,pormap,bcflag,courant,totaltime,kmap,wells,porousarea,ds,tconv,...
        fract,elem_post,regionfract,coord_post,formethod,timedimkey,phflw,knownb,...
        rowposit,region,numsline,substeps,maxsteps,jointnodes] = preprocessor(inFile)
    
global afrat kfrat

%--------------------------------------------------------------------------
%It reads the problem data file
pdata = INI('File', inFile); pdata.read();

% -------
% problem
filepath = char(pdata.problem.meshFilePath);  %Define the file path
res_folder = char(pdata.problem.outputPath);
if exist(res_folder,'dir') == 0
    mkdir(res_folder);
end

formethod = pdata.problem.pSolver;
transolver= pdata.problem.tSolver;
if strcmp(transolver,'NONE')
    phflw=1; % only pressure problem to solve
else
    phflw=2; % pressure and transport to solve
    if strcmp(transolver,'IMPLICIT')
        ds=1;
    elseif strcmp(transolver,'EXPLICIT')
        ds=2;
    elseif strcmp(transolver,'STREAMLINE')
        ds=3;
    end
end


% -----------------
% fluids properties
dens = pdata.fluids.density;
visc = pdata.fluids.vicosity;
satlimit = pdata.fluids.residualSaturation;

% --------------
% rock propertie
buf = pdata.rock.permeability; [ibuf, ~] = size(buf);
kmap = zeros(ibuf, 5); kmap(:,2:end) = buf; kmap(:,1) = 1:ibuf;
pormap = pdata.rock.porosity;

% -------------------
% boundary conditions
bcflag = pdata.boundaryCondition.bcflagAtFaces;

% ----------------------------------------
% Convergence criteria and simulation time
courant = pdata.transportSolver.courant;
timedimkey = 's';
totaltime = [0; pdata.problem.totalTime];

% -----
% wells
buf = pdata.source.injWellLocation;
inj = []; numwell = 0;
if ~isnan(buf)
	[ibuf, ~] = size(buf);
	assert(ibuf == numel(pdata.source.injWellRadius), ...
           'incompatible sizes between injWellLocation and injWellradius');
	assert(ibuf == numel(pdata.source.injWellsat), ...
           'incompatible sizes between injWellLocation and injWellsat');
    sat = pdata.source.injWellsat;
	assert(ibuf == numel(pdata.source.injWellPressure), ...
           'incompatible sizes between injWellLocation and injWellPressure');
	assert(ibuf == numel(pdata.source.injWellRate), ...
           'incompatible sizes between injWellLocation and injWellRate');
    jflag = 401*ones(ibuf,1); jflag(isnan(pdata.source.injWellPressure)) = 0;
    value = pdata.source.injWellPressure;
    if any(isnan(pdata.source.injWellPressure))
        value(isnan(pdata.source.injWellPressure)) = ...
                pdata.source.injWellRate(isnan(pdata.source.injWellPressure));
    end
	assert(ibuf == numel(pdata.source.injWellType), ...
           'incompatible sizes between injWellLocation and injWellType');
    wtype = pdata.source.injWellType;
    iflag = 301*ones(ibuf,1);
    
    %      id_injection, x y, well influence radius, 301 flags, inj sat, pressure flag, pressure flag value, well type
    inj = [ones(ibuf,1), buf, pdata.source.injWellRadius', iflag, sat', jflag, value', wtype'];
    numwell = ibuf;
end


buf = pdata.source.prdWellLocation;
prd = [];
if ~isnan(buf)
	[ibuf, ~] = size(buf);
	assert(ibuf == numel(pdata.source.prdWellRadius), ...
           'incompatible sizes between prdWellLocation and prdWellradius');
	assert(ibuf == numel(pdata.source.prdWellPressure), ...
           'incompatible sizes between prdWellLocation and prdWellPressure');
	assert(ibuf == numel(pdata.source.prdWellRate), ...
           'incompatible sizes between prdWellLocation and prdWellRate');
    jflag = 501*ones(ibuf,1); jflag(isnan(pdata.source.prdWellPressure)) = 0;
    value = pdata.source.prdWellPressure;
    if any(isnan(pdata.source.prdWellPressure))
        value(isnan(pdata.source.prdWellPressure)) = ...
                pdata.source.prdWellRate(isnan(pdata.source.prdWellPressure));
    end
	assert(ibuf == numel(pdata.source.prdWellType), ...
           'incompatible sizes between injWellLocation and injWellType');
    wtype = pdata.source.prdWellType;
    iflag = zeros(ibuf,1);
    sat = zeros(ibuf,1);
    
    %      id_injection, x y, well influence radius, 301 flags, prd sat, pressure flag, pressure flag value, well type
    prd = [2*ones(ibuf,1), buf, pdata.source.prdWellRadius', iflag, sat, jflag, value', wtype'];
    numwell = numwell + ibuf;
end
well = [inj; prd];


%-------------------------------------------------------------------------%
%-------------------------- Fracture Options -----------------------------%
afrat = 0; kfrat=0; ntypefrat=0; hybridgrid=0; numsline = 0; substeps = 0;
maxsteps = 0; tconv = 0;
if isprop(pdata,'fracture')
    hybridgrid=1;
    ntypefrat = pdata.fracture.nfrac;
	afrat = pdata.fracture.thickness;
	kfrat = pdata.fracture.kfracture;
end

if and(isprop(pdata,'streamline'), ds==3)
    numsline = pdata.streamline.numsline;
    substeps = pdata.streamline.substeps;
    maxsteps = pdata.streamline.maxsteps;
end
%-------------------------------------------------------------------------%

%Jump a line in the matlab's prompt
disp(' ');
%Call atention to data to be generated
disp('---------------------------------------------------');
disp('>> It is generating the Data Structure...');
disp(' ');


%--------------------------------------------------------------------------
%"coord" matrix - The coordinates (x,y,z) of each point in a carts. 

[coord,nnode] = getcoord(filepath);

%Gives the information of "coord" generated
disp('"coord" was generated!');

%--------------------------------------------------------------------------
%"elem" matrix - The points that constitute each element

[elem,nbe,nelem,nodelim,intnode,flaglim,fract,nodefrat,...
      spclnodes,region,regionfract] = getelem(filepath,nnode,numwell,well);

%Fill the matrix "centelemcoord" ------------------------------------------
centelem = getcentelem(coord,elem);
[ elem ] = reordnodelem( elem, centelem, coord );
%--------------------------------------------------------------------------

elem_post=elem;
coord_post=coord;
jointnodes=0;

% Parcial do bedge -------------------------------------------------------%
%Open the *.msh file
readbound = fopen(filepath);
%In principle "bedgeaux" read in its fourth column the geometry flag from 
%*.msh file generated by gmsh. After that "bedge" receives this info.
getboundata = textscan(readbound,'%*u %*u %*u %u %*u %u %u',nbe,...
    'HeaderLines',8 + nnode + nodelim + nodefrat + intnode);
%Attribute the data to "bedgeaux"
bedgeaux(:,5) = getboundata{1};
bedgeaux(:,1:2) = [getboundata{2} getboundata{3}];
%Close file
fclose(readbound);
%-------------------------------------------------------------------------%

%Hybrid-grid construction
if hybridgrid==1 && ntypefrat>0
    
    if (sum(elem(:,4)==0)>0)&&(ds==3)
       ds=1;
       disp('---------------------------------------------------');
       disp('N?o foi poss?vel utilizar o m?todo das streamlines!');
       disp('Seguiremos com o SEQ cl?ssico!');
       disp('Para streamlines, utilize malhas quadrilaterais!');
       disp('---------------------------------------------------');
    end
    
%     if ds==3
%           
          [ coord, elem, fract, kmap, region, regionfract, flagcorr ] = quadhybgrid_build( coord, elem, ...
                         centelem, fract, kfrat, kmap, spclnodes, region, regionfract, flaglim );
    [~,j] = size(fract);
    if j > 4
        fract = fract(:,1:4);
    end
%                                  
%     else
%                                  
%         [ coord, elem, fract, kmap, region, regionfract ] = hybridgrid_build( coord, elem, ...
%                fract, afrat, kfrat, kmap, spclnodes, region, regionfract );
%                      
%           [ coord, elem, fract, kmap, region, regionfract ] = tuliogrid_build( coord, elem, ...
%                        fract, kfrat, kmap, spclnodes, region, regionfract, flaglim );
%           flagcorr = 0;
%         
%     end
       
    elem(nelem+1:nelem+size(fract,1),1:size(fract,2))=fract;
    region(nelem+1:nelem+size(fract,1))=regionfract;
    fracttype = unique(regionfract);
    pormap(size(pormap,1)+1:size(fracttype,1)+size(pormap,1),1)=ones(size(fracttype,1),1);
    nelem=size(elem,1);
    nnode=size(coord,1);
    jointnodes=spclnodes(spclnodes(:,2)>2000,3);
    
else
    flagcorr = 0;
end
vtkWriter([res_folder,filesep,'meshLDFM'],'coord',coord, 'elem', elem,...
    'timestep', 0, 'cycle', 0, 'Regions', elem(:,end));

%Gives the information of "elem" generated
disp('"elem" was generated!');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%"centelem" matrix - It get the centroid coordinate for each control volume 

%Fill the matrix "centelemcoord"
centelem = getcentelem(coord,elem);

%Gives the information of "elem" generated
disp('"centelem" was generated!');

%--------------------------------------------------------------------------
[ elem ] = reordnodelem( elem, centelem, coord );
%--------------------------------------------------------------------------
%"elemarea" vector - It get the area for each control volume 

%Fill the vector "elemarea"
[elemarea] = calcelemarea(coord,elem);

%Gives the information of "elem" generated
disp('"elemarea" was generated!');

%--------------------------------------------------------------------------
%"bedge" and "inedge" - It get inform. about boundary and internal edges

[bedge,inedge,knownb] = getinedgebedge(coord,elem,bedgeaux,flaglim,flagcorr);

[inedge,bedge] = reordbedge(inedge,bedge,coord,centelem);
                    
%Gives the information of "bedge" and "inedge" generated
disp('"bedge" was generated!');
disp('"inedge" was generated!');

%"normals" is a matrix with the normal vectors (normal to left element)

%Get the normals
normals = calcnormals(coord,centelem,bedge,inedge);

%Gives the information of "normals" generated
disp('"normals" was generated!');

%--------------------------------------------------------------------------
%"esurn" vectors - Report the amount of ELEMENTS surronding each POINT.

%This is constituted by two vectors as a linked list. The first one ...
%("esurn1") contain the number of each element associated with the node ...
%listeds. The second one ("esurn2") contain numbers related to position...
%(on "esurn1" vector) of start of counter.

%Swept the "Elem" matriz from 1 until Nelem (rows) and inside it swept from
%1 until the number of nodes that constitute each element (columns)

%Create and initialize "esurn2"
%nnode + 1 rows and 1 column. Begins from 0 until "nnode"
esurn2 = zeros(nnode + 1,1); 
%Ratify "elem" size
nelem = size(elem,1);

    %"ielem" is a counter of elements (swept the rows of "elem")
    for ielem = 1:nelem
        %Get "elemcontents"
        elemcontents = elem(ielem,1:5);
        %Get the amount of edges by element.
        amountedges = sum(elemcontents ~= 0);
        %Get a vector with the vertices for the element evaluated. 
        vertices = elemcontents(1:amountedges);
        %Get the number of vertices for this element.
        numvert = length(vertices);
        %This loop swept all nodes who constitute each element (columns of 
        %"elem" matrix). 
        for inode = 1:numvert    
            %"nodeval" is the receives the number of each node. 
            %The elem(1,1) has nodeval = 1, elem(1,2) has nodeval = 5 
            %(see basic exemple in Marcio's notebook) 
            nodeval = vertices(inode);
            %Add one to position which corresponds to number of node. We
            %have to consider that "esurn2" start with value 0 for a better
            %use as will be saw in the future (see Lohner, 2001 - Chap 2)
            esurn2(nodeval + 1) = esurn2(nodeval + 1) + 1;
        end  %End of FOR (in each element)
    end  %End of FOR (in all "elem" matrix)

%"esurnqnt" recives "esurn2" before this be recalculated in order to 
%represent an acumulate value (function "cumsum") 
esurnqnt = esurn2;    
%Create a vector with the amount accumulated from "esurn2" already 
%obtained (esurn2cum(i) = esurn2(i) + esurn2(i-1))
esurn2 = cumsum(esurn2);
esurn2aux = esurn2;

%Create and initialize "esurn1"
%The number of rows of "esurn1" is equal to amount accumulated of elements
%in each node
esurn1 = zeros(esurn2(nnode + 1),1);
%Initialize "nodestore" which control if the node evaluated 
%already was verifyed. If the node in the "elem" was already evaluated 
%inode is added in one.
nodestore = 0;
%"inodestore" is a counter of "nodestore"
inodestore = 1;

%Fill "esurn1" (list of elements surrounding each node)
    %Swept from 1 until the last element
    for ielem = 1:nelem  
        %Get "elemcontents"
        elemcontents = elem(ielem,1:5);
        %Get the amount of edges by element.
        amountedges = sum(elemcontents ~= 0);
        %Get a vector with the vertices for the element evaluated. 
        vertices = elemcontents(1:amountedges);
        %Get the number of vertices for this element.
        numvert = length(vertices);
        %Swept each node of row's "elem" evaluated
        for inode = 1:numvert   
            %"nodeval" receives the index of node on "elem" matrix
            nodeval = vertices(inode);
            %Just one element concors to node evaluated.
            %In general this occur in the domain corner.
            if esurnqnt(nodeval + 1) == 1
                %"store" obtain the value (of "esurn2aux" vector) which 
                %corresponds to position whose number is "nodeval" + 1.
                %Remenbering "esurn2" begins from 0 (for that +1)
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "Elem" row in the "esurn1"
                %whose the position is termined by "strore" index
                esurn1(store) = ielem;
                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - 1;
            %Two or more elements concor to node evaluated and the node 
            %evaluated to be inside domain. In general is what happen in 
            %the domain.
            elseif esurnqnt(nodeval + 1) > 1 && ...
                    ismember(nodeval,nodestore) ~= 1
                %Initialize "elemorder". This parameter receives the
                %elements which try to node evaluated in increasing order.
                elemorder = zeros(1,esurnqnt(nodeval + 1));
                
                %----------------------------------------------------------
                %"nodeval" over boundary
                
                %Verify if the node evaluated ("nodeval") is over the 
                %boundary. In this case we have to find the first element
                %and the last element in order maintain the counter
                %clockwise way (little more complicated).
                rowsinbedge = ...
                    bedge(any(logical(bedge(:,1:2) == nodeval),2),1:3);
                %The node evaluated is inside domain. In this case is
                %easer than the earlier situation once any element may be 
                %the first one. The cycle finishes when the last element 
                %before the first is found.
                if isempty(rowsinbedge)  
                    %Attribute the first element to "elemorder". 
                    %Considering "ielem" contain one of elements evaluated.
                    elemorder(1) = ielem;

                %The node evaluated is over the boundary. The first element
                %in "elemorder" cames from "bedge"
                else
                    %Initialize "vertexout" and "firstelem"
                    vertexout = 0;
                    firstelem = 0;
                    %Initialize "elemcandidate". It stores the number of
                    %elements which can be the first element in "elemorder"
                    elemcandidate = rowsinbedge(:,3);
                    %"rowsvertices" receives just the number of vertices
                    rowsvertices = rowsinbedge(:,1:2);
                    %We need verify which other vertices belong to edge
                    %over the boundary.
                    vertoutnodeval = setdiff(reshape(rowsvertices',4,1),...
                        nodeval,'stable');
                    
                    %Evaluate which edge is the first:
                    %Get the vector (first edge)
                    edgevec = coord(vertoutnodeval,:) - ...
                        [coord(nodeval,:); coord(nodeval,:)];
                    %Get a vector which points to out of domain.
                    pointout = [coord(nodeval,:); coord(nodeval,:)] - ...
                        centelem(elemcandidate,:);
                    
                    %Attribute the other vertex (edge over the boundary).
                    vertexout = vertexout + vertoutnodeval(1)*...
                        (cross(pointout(1,:),edgevec(1,:)) > 0);
                    vertexout = vertexout + vertoutnodeval(2)*...
                        (cross(pointout(2,:),edgevec(2,:)) > 0);
                    %Attribute the first element sharing the edge over the 
                    %boundary.
                    firstelem = firstelem + elemcandidate(1)*...
                        (setdiff(cross(pointout(1,:),edgevec(1,:)),0) > 0);
                    firstelem = firstelem + elemcandidate(2)*...
                        (setdiff(cross(pointout(2,:),edgevec(2,:)),0) > 0);

                    %Attribute the first element to "elemorder(1)"
                    elemorder(1) = firstelem;
                end  %End of internal IF

                %Create a vector with number of "inedge" rows
                iinrow = 1:size(inedge,1);

                %Find the last element to be stored in the "elemorder"
                for iorder = 2:esurnqnt(nodeval + 1)
                    %Initialize "rowpointer". This parameter must be
                    %initialized in order to avoid superposition of values.
                    rowpointer = 0;
                    %Obtain the row or rows of "inedge" whose "nodeval" 
                    %belongs to first two columns of "inedge" and, to same 
                    %time, the aforefind element belongs to two last 
                    %columns of "inedge".
                    %Evaluate the vertices and elements in "inedge":
                    evalvertelem = all([any(ismember(inedge(:,1:2),...
                        nodeval),2) any(ismember(inedge(:,3:4),...
                        elemorder(iorder-1)),2)],2);
                    %Get the "inedge" row.
                    rowinedge = iinrow(evalvertelem)';

                    %This loop points to entity of "rowinedge" which
                    %satisfy the necessary conditions. That row which
                    %satisfies receive the value "1", while that one 
                    %which does not satisfy receive "0"
                    for irow = 1:length(rowinedge)
                        rowpointer(irow) = ...
                            (inedge(rowinedge(irow),2) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),3)) | ...
                            (inedge(rowinedge(irow),1) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),4));
                    end  %End of internal FOR

                    %"rowdef" is a definitive row which satisfy the
                    %necessary conditions aforesaid. That is the 
                    %position which have "1" insteady "0" 
                    rowdef = logical(rowpointer == 1);   
                    %"rowinedge" receives, finaly, the unic row which 
                    %satisfy the conditions above.
                    rowinedge = rowinedge(rowdef);
                    %A decision is made as function of "nodeval"'s 
                    %value. The commands written below avoids "IF"  
                    %"nodeval" is lower than other node which 
                    %constitute the common edge among two elements. 
                    %In this case, a left element is attributed to 
                    %"elemorder" just if "nodeval" is lower than other 
                    %element
                    rowinedge=min(rowinedge);
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,3)*any((nodeval ~= ...
                        max(inedge(rowinedge,1:2))));
                    %In this case, a right element is attributed to 
                    %"elemorder" just if "nodeval" is major than other 
                    %element   
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,4)*any((nodeval == ...
                        max(inedge(rowinedge,1:2))));
                end  %End of FOR
                
                %"nodestore" receives the node evaluated. Thus, when this
                %will evaluated in the future that will be by passed,
                %because was already (see line 1241)
                nodestore(inodestore) = nodeval;
                %Add 1 to "inodestore"
                inodestore = inodestore + 1;

                %Put the number in its places in the "esurn1"
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "elem" row in the 
                %"esurn1" whose the position is termined by "strore" index
                
                %"retrocount" is a retroative counter.
                retrocount = esurnqnt(nodeval + 1);
                for iretro = store:-1:store - (esurnqnt(nodeval+1) - 1)
                    esurn1(iretro) = elemorder(retrocount);
                    %backward the counter
                    retrocount = retrocount - 1;
                end  %End of FOR

                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - ...
                    esurnqnt(nodeval + 1);
            end  %End of external IF
            %Jump a node position if any conditions above is atended
            inode = inode + 1;  %Next node
        end  %End of WHILE (read each column of each row of "elem")
    end  %End of FOR (read each row of "elem")
    
%Gives the information of "esurn" generated
disp('"esurn" was generated!');

%--------------------------------------------------------------------------
%"nsurn" vectors - Report the amount of NODES surronding each NODE.

%We will use the vectors "eserp1" and "esurn2" to swept the elements around
%each point. In the evaluated element we must find the near points to point
%who corresponds to position of "Esurn2" (the first position of this vector 
%corresponds to point 1 and consecutively)

%Initialization of variables
%Works as a counter to "nsurn1" vector
% store = 0;
%"nsurn2" points to "nsurn1" and denote the begin and end of node 
%surrounding cycle of each node

nsurn2 = zeros(nnode + 1,1);

%Swept the vector "esurn2" to pass for each point (first for). Each one 
%index "i" is a position of aforesaid (abovementioned) vector.
%With the development of "i", points already visited will not any more  
    for i = 1:nnode
        %Swept from position "i" add 1 until position "i+1" (end of 
        %elements surrounding point cycle).
        %"iesurn" is a counter of thous elements surrounding each point
        store = 1;
        %Initialization of variable
        %"nsurn1aux" is an auxiliar vector which avoid repeated node 
        %be stored in "nsurn1" (see below)
        nsurn1aux = 0;
        %"insurn" is a counter who will be used in "nsurn1aux" vector
%         insurn = 0;
        %Swept all elements surround each node evaluated
        for iesurn = esurn2(i) + 1:esurn2(i + 1)
            %"ielem" will count which elements be in the cycle.
            %This parameter receives the elements 1,2 & 3 (whose positions 
            %are 7,8 and 9) - see basic example in Marcio's notebook
            
            ielem = esurn1(iesurn);   
            
            numn = sum(elem(ielem,:)~=0);
            posn = find(elem(ielem,:)==i);
                       
            refn = zeros(1,numn+2);            
            refn(2:numn+1) = 1:numn;
            posr = find(refn==posn);
            
            refn(numn+2) = 1; refn(1) = numn; 
            
            vecn = [elem(ielem,refn(posr-1)) elem(ielem,refn(posr+1))];
            
            for j=1:size(nsurn1aux,2)
                if nsurn1aux(j)==vecn(1)
                    vecn(1)=0;
                end
                if nsurn1aux(j)==vecn(2)
                    vecn(2)=0;
                end
            end
            rsn = sum(vecn~=0);  
            
            if rsn==2
                vsn1 = coord(vecn(1),:)-coord(i,:);
                vsn2 = coord(vecn(2),:)-coord(i,:);
                dsn = cross(vsn1,vsn2);
                if dsn(3)<0
                    dsn(1)=vecn(2);vecn(2)=vecn(1);vecn(1)=dsn(1);
                end
            end
            
            nsurn1aux(store:store+rsn-1)=nonzeros(vecn);
            store = size(nsurn1aux,2)+1;
  
        end  %End of second FOR
        %Update "nsurn2". This variable receives the number of increments
        %that "store" had
        nsurn2(i+1) = nsurn2(i) + store-1;
        nsurn1(nsurn2(i)+1:nsurn2(i+1),1) = nsurn1aux';
    end  %End of first FOR

%Reorder "nsurn1" (counterclockwise). It returns a column vector.
% [nsurn1,nsurn2] = reordernsurn(esurn1,esurn2,nsurn1,nsurn2,bedge,...
%                     coord,elem);
                
rowposit = getrowposition(bedge,inedge,nsurn1,nsurn2);

%Gives the information of "nsurn" generated
disp('"nsurn" was generated!');
%-------------------------------------------------------------------------%
        
%--------------------------------------------------------------------------
%Generate "esureface" (face neighbor element) and "esurefull" (vertices and 
%face) neighbor

[esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,inedge,...
    esurn1,esurn2,nsurn1,nsurn2);

%Gives the information of "nsurn" generated
disp('"esure" were generated!');


disp('physical properties were generated!');

for j=1:size(elem,1)
   porousarea(j)=elemarea(j)*pormap(region(j)); 
end
kfrat=kfrat(:,1);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%          WELLS - Gives details about well(s) set on the domain          %
%-------------------------------------------------------------------------%
%Look to line 372. There the data referred to wells (written in the 
%"Start.dat") are read and used in this section.

%Initialize "wells" and "countwell". 
wells = 0;
countwell = 0;
%"j" is a conter of wells row.
j = 0;
%When the WELL is a POINT (Unstructured mesh with mesh adaptation)
%Find the nodes which contain the wells. These are related with rows of
%matrix "coord".
if numwell > 0 && any(well(:,9) == 0)
    %Swept the amount of wells
    for iwell = 1:numwell
        %"centwell" is a vector with the node which stores the center 
        %of well set. Surrounding these nodes the elements will have 
        %sourse terms.
        centwell(iwell) = find(coord(:,1) == well(iwell,2) & ...
            coord(:,2) == well(iwell,3));
        %Once the node was finded we can verify which elements to be
        %surrounding it. This information will be stored in the first
        %column of "wells matrix"
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,1) = esurn1(esurn2(centwell(iwell))+1:...
            esurn2(centwell(iwell)+1));
        %The second column will store the well's number. This order
        %follow the order of write in the "Start.dat". "iwell" gives
        %this number.
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,2) = iwell;
        %The third column stores the saturation flag. This may have 
        %values from [301 - 400]. Is important to remember that to 
        %productor well this parameter receives "0".             
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,3) = well(iwell,5);
        %The fourth column stores the saturation value. This will be "0"
        %for producer wells.
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,4) = well(iwell,6);
        %The fiveth column stores the pressure flag. That may be have 
        %values among [401 - 500] to injector well and [501 - 600] to 
        %productor well  
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,5) = well(iwell,7);
        %The sixth column stores the pressure value. That may receives the
        %flow rate value if the pressure flag is "0"
        wells(1 + countwell:...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell))) + ...
            countwell,6) = well(iwell,8);
        %"countwell" receives a new value
        countwell = countwell + ...
            (esurn2(centwell(iwell)+1) - esurn2(centwell(iwell)));
    end  %End of FOR
end  %End of IF

%When the WELL is a POINT (STRUCTURED mesh with NO mesh adaptation)
if numwell > 0 && any(well(:,9) == 1)
    %Verify which ones row of "well" has sign related to point well in 
    %structured mesh (1)
    numstrwell = find(well(:,9) == 1);
    %Cont the number of wells
    countwell = countwell + 1;
    %Swept the wells alocated.
    for istrwell = 1:length(numstrwell)
        %Verify in "coord" the coordinate "x" of well center
        fstnodex = logical(coord(:,1) == well(numstrwell(istrwell),2)); 
        %Verify in "coord" the coordinate "y" of well center
        fstnodey = logical(coord(:,2) == well(numstrwell(istrwell),3)); 
        %Obtain the row of "coord" where the coordinate of "x" and "y" are
        %coincident. It is the number of node where the well is centered.
        wellnode = find((fstnodex + fstnodey) == 2);
        %Verify the elements which surrounding this node.
        %"poselem" guards the position in "esurn1" of elements which
        %surrounding the node  "wellnode"
        poselem = esurn2(wellnode) + 1:esurn2(wellnode + 1);

        %Fill the first column of "wells" (elements)
        wells(j + 1:j + length(poselem),1) = esurn1(poselem);
        %Fill the second column of "wells" (number of well)
        wells(j + 1:j + length(poselem),2) = countwell;
        %Fill the thrid column of "wells" (flag of saturation)
        wells(j + 1:j + length(poselem),3) = well(numstrwell(istrwell),5);
        %Fill the fourth column of "wells" (value of saturation)
        wells(j + 1:j + length(poselem),4) = well(numstrwell(istrwell),6);
        %Fill the fiveth column of "wells" (flag of pressure)
        wells(j + 1:j + length(poselem),5) = well(numstrwell(istrwell),7);
        %Fill the sixth column of "wells" (value of pressure)
        wells(j + 1:j + length(poselem),6) = well(numstrwell(istrwell),8);
        %Fill the thrid column of "wells" (coord 1)
        wells(j + 1:j + length(poselem),7) = well(numstrwell(istrwell),2);
        %Fill the fourth column of "wells" (coord 2)
        wells(j + 1:j + length(poselem),8) = well(numstrwell(istrwell),3);
        %Increment "j"
        j = j + length(poselem);
        %Increment "countwell" (when there is other well)
        countwell = countwell + 1;
    end  %End of FOR        
end  %End of IF

%BUCKLEY-LEVERETT Applications.
%When the WELL is a LINE and this LINE is over BOUNDARY
if numwell > 0 && any(well(:,9) == 2)
	%Verify which ones row of "well" has sign related to BL well (2)
	blwell = find(well(:,9) == 2);
	%Catch the elements associated with edge defined.
    xbfaces = [mean(reshape(coord(bedge(:,1:2),1), size(bedge,1),2),2),...
               mean(reshape(coord(bedge(:,1:2),2), size(bedge,1),2),2)];
	%This "for" conts each two points due to each line is defined by two
	%points
	for ibl = 1:2:length(blwell)
		%Cont the number of wells
		countwell = countwell + 1;
		
        % line on boundary
        xline = well(blwell(ibl),2:3) - well(blwell(ibl + 1),2:3);
        idx = find(xline == 0);
        brow = xbfaces(:,idx) == well(blwell(ibl),1+idx);
        
        %Fill the first column of "wells" (elements)
        wells(j + 1:j + sum(brow),1) = bedge(brow,3);
        %Fill the second column of "wells" (number of well)
        wells(j + 1:j + sum(brow),2) = countwell;
        %Fill the thrid column of "wells" (flag of saturation)
        wells(j + 1:j + sum(brow),3) = well(blwell(ibl),5);
        %Fill the fourth column of "wells" (value of saturation)
        wells(j + 1:j + sum(brow),4) = well(blwell(ibl),6);
        %Fill the fiveth column of "wells" (flag of pressure)
        wells(j + 1:j + sum(brow),5) = well(blwell(ibl),7);
        %Fill the sixth column of "wells" (value of pressure)
        wells(j + 1:j + sum(brow),6) = well(blwell(ibl),8);
        %Increment "j"
        j = j + sum(brow); 
    end  %End of FOR
end  %End of IF

%When the WELL is a LINE 
%There is a source line INSIDE domain
if numwell > 0 && any(well(:,9) == 3) 
    %Define the well number
    wellnumber = countwell + 1;
    for i = 1:size(inboundedge,1)
        %Find the row of "inedge" which contain the two nodes which
        %constitute the line source
        indexedge = logical(inedge(:,1) == inboundedge(i,1) & ...
            inedge(:,2) == inboundedge(i,2));
        %Once the edge was finded we can verify which elements to be
        %surrounding it. This information will be stored in the first
        %column of "wells matrix"
        %Attribute the lft and right elements
        wells(countwell + i:countwell + i + 1,1) = ...
            [inedge(indexedge,3); inedge(indexedge,4)];
        %The second column will store the well's number. This order
        %follow the order of write in the "Start.dat". "iwell" gives
        %this number.
        wells(countwell + i:countwell + i + 1,2) = wellnumber;
        %The third column stores the saturation flag. This may have 
        %values from [301 - 400].              
        %Find a row of "well" which contain 1 in the seventh column.
        wellrow = logical(well(:,9) == 3);
        %Fill the third column of "wells" (flag of saturation)
        wells(countwell + i:countwell + i + 1,3) = well(wellrow(1),5);
        %Fill the fourth column of "wells" (value of saturation)
        wells(countwell + i:countwell + i + 1,4) = well(wellrow(1),6);
        %The fiveth column stores the pressure flag. That may be have 
        %values among [401 - 500] to injector well and [501 - 600] to 
        %productor well  
        wells(countwell + i:countwell + i + 1,5) = well(wellrow(1),7);
        %Fill the sixth column of "wells" (value of pressure)
        wells(countwell + i:countwell + i + 1,6) = well(wellrow(1),8);
        %Increment "countwell"
        countwell = countwell + 1;
    end  %End of FOR
end  %End of IF

%The domain has no wells
if numwell == 0
    wells = 0;
end  %End of IF

%Message of preprocessor.
disp('wells properties were generated!');
%Jump a line
disp(' ');
%Final message
disp('>> All data were generated with success!!!');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%FUNCTION "getcoord"
%--------------------------------------------------------------------------
function [coord,nnode] = getcoord(filepath)
%Open the *.msh file
readmsh = fopen(filepath);

%"nnode" is the number of nodes in the discrete domain
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',4);
%Attribute the data to "nnode"
nnode = getmshdata{1};

%"coord" is a matrix which contain the coordinate of each point of domain
getmshdata = textscan(readmsh,'%*u %f64 %f64 %f64',nnode,'HeaderLines',1);
%Fill the matrix "coord" with the x, y and z coordinates.
coord = cell2mat(getmshdata);

%Close the *.msh file
fclose(readmsh);

%--------------------------------------------------------------------------
%FUNCTION "getelem"
%--------------------------------------------------------------------------
function [elem,nbe,nelem,nodelim,intnode,flaglim,fract,nodefrat,...
         spclnodes,region,regionfract] = getelem(filepath,nnode,numwell,well)
%Identify on the *.msh file a like-element type. Actualy, the *.msh
%gives entities which hold points, edges and elements (like-elements).
%This entities will all of them be reads.
%Open the *.msh file
readmsh = fopen(filepath);
        
%"nent" is the number of all entities understud as element in *.msh file
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
%Attribute the data to "nnode"
nent = getmshdata{1};
        
%Once we have "nent", both like-element and external edges entities may be 
%located. "meshtype" is the type of element used in domain discretization.
%2 ==> triangles; 3 ==> quadrangles. But that will be used also to define
%how many points entities there on the *.msh file

%There is an INTERNAL LINE
if numwell > 0  && ( any(well(:,9) == 3) || any(well(:,9) == 0) )
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n%n%*n%n%*n%n%n%*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Gives the contribution to "entitype"
    entitype = auxmat(:,1);
    entiflag = auxmat(:,2);
    %Get "inboundedge" from "auxmat"
    inboundedge = ...
        auxmat(logical(auxmat(:,1) ~= 15 & auxmat(:,2) > 1000),3:4);
    %Orient the order in each column
    inboundedge = sort(inboundedge,2);
    
    %"nbe" verifies how many entities are edges which constitute the 
    %boundary
    nbeaux = sum(entitype == 1)-sum(entiflag > 2000);
    nbe = nbeaux - size(inboundedge,1);
    nfrat = sum(entiflag > 2000);
    %"nodelim". These nodes constitute the limits of geometry. 
    %Get the amount of internal node
    intnode = sum(logical(auxmat(:,1) == 15 & auxmat(:,2) > 1000));
    %We exclude any node inside the domain.    
    nodelim = sum(entitype == 15) - intnode;
%Other cases
else
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Fill "entitype"
    entitype = auxmat(:,1);
    entiflag = auxmat(:,2);
    %"nbe" verifies how many entities are edges which constitute the 
    %boundary
    nbe=0;
    nbeaux=0;
    nfrat=0;
    for ie=1:size(entitype,1)
        if entitype(ie)==1 && entiflag(ie)<2000 && entiflag(ie)>100
            nbe=nbe+1;
        elseif entitype(ie)~=15 && entiflag(ie)>2000
            nfrat=nfrat+1;
        end
        if entitype(ie)==1 && entiflag(ie)<2000
            nbeaux=nbeaux+1;
        end
    end
    %"nodelimit". These nodes constitute the limits of geometry
    nodelim = sum(entitype == 15);
    intnode = 0;
end  %End of IF

%Define the amount of each entity:
%Verify how many "entities" are nodes. This value is attributed to 
%"ntri" verifies how many triangles there are in the domain
ntri = sum(entitype == 2); 
%"nquad" verifies how many quadrangles there are in the domain
nquad = sum(entitype == 3);
%"nelem" is the number of elements in the domain (triangles and/or 
%quadrangles)
nelem = ntri + nquad + nfrat;

%Get "flaglim" from "auxmat". It is a vector with the flag associated to
%external nodes.
flaglim = auxmat(logical(auxmat(:,1) == 15 & auxmat(:,2) < 1000),2);

%Close the *.msh file
fclose(readmsh);

%---------------------------
%Define the matrix ("elem"):

%Open again the *.msh file
readmsh = fopen(filepath);

%"elem" is a matrix which contain the nodes which constitute it
%Create and initialize the parameter "elem"
elem = zeros(nelem,5);  %nelem rows and 5 columns
region = zeros(nelem,1);  %region flag

%Fill the matrix "elem" with 4 nodes for a quadrangle element and 3 nodes 
%for a triangle element. In this case the fourth column receives null
%value. The fiveth column receive material propertie flag.
if nfrat > 0  
    %Fractures
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u',nfrat,...
        'HeaderLines',8 + nnode + nodelim + intnode + nbeaux);
    elemaux = cell2mat(getmshdata);
    %Region flag
    region(1:nfrat) = elemaux(:,1);
    %Attribute to "elem" (another column)
    elem(1:nfrat,1:2) = elemaux(:,2:3);
end  %End of first IF
if ntri > 0  
    %Triangle
    if nfrat==0
        H1=8 + nnode + nodelim + intnode + nbeaux;
    else
        H1=0;
    end
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u',ntri,...
        'HeaderLines',H1);
    elemaux = cell2mat(getmshdata);
    %Region flag
    region(nfrat+1:nfrat+ntri) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(nfrat+1:nfrat+ntri,1:3) = elemaux(:,2:4);
end  %End of first IF
if nquad > 0  
    %Quadrangles
    if (nfrat==0)&&(ntri==0)
        H1=8 + nnode + nodelim + intnode + nbeaux;
    else
        H1=0;
    end
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u %u',nquad,...
        'HeaderLines',H1);
    elemaux = cell2mat(getmshdata);
    %Region flag
    region(ntri+nfrat+1:ntri+nfrat+nquad) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(ntri+nfrat+1:ntri+nfrat+nquad,1:4) = elemaux(:,2:5); 
end  %End of second IF

fract=elem(1:nfrat,:);
regionfract=region(1:nfrat);
elem1=elem(nfrat+1:end,:);
region1=region(nfrat+1:end);
clear elem region
elem=elem1;
region=region1;
nelem=nelem-nfrat;

nlim=0;
nodefrat=0;
for in=1:nodelim
    if entitype(in)==15 && entiflag(in)<2000
       nlim=nlim+1; 
    elseif entitype(in)==15 && entiflag(in)>2000
       nodefrat=nodefrat+1; 
    end
end

nodelim=nlim;
%Close the *.msh file
fclose(readmsh);

% Open it again to get information about the fracture nodes.
readmsh = fopen(filepath);      
getmshdata = textscan(readmsh,'%*u %u %*u %u %*u %u',nodelim + nodefrat,'HeaderLines',8 + nnode);
spclnodes = cell2mat(getmshdata);
fclose(readmsh);

%--------------------------------------------------------------------------
%Function "calcnormals"
%--------------------------------------------------------------------------
function [normals] = calcnormals(coord,centelem,bedge,inedge)
%Define the matrix rotation
R = zeros(3);
R(1,2) = 1;
R(2,1) = -1;

%Initialize "normals"
normals = zeros(size(bedge,1) + size(inedge,1),3);

%Fill the matrix "normals" (referred to edges in boundary)
for iflowb = 1:size(bedge,1)
    %Fill "normals"
    normals(iflowb,:) = ...
        R*(coord(bedge(iflowb,2),:) - coord(bedge(iflowb,1),:))';

    %Confirm the normal orientation (out of control volume)
    %Calculate a pointer (centroid to out of control volume)
    pointer = coord(bedge(iflowb,1),:) - centelem(bedge(iflowb,3),:);
    %Correct (if necessary) the normal orientation
    normals(iflowb,:) = ...
        normals(iflowb,:)*sign(dot(pointer,normals(iflowb,:)));
end  %End of FOR ("bedge")

%Fill matriz "normals" (referred to edges into domain)
for iflowin = 1:size(inedge,1)
    %Fill "normals"
    normals(size(bedge,1) + iflowin,:) = ...
        R*(coord(inedge(iflowin,2),:) - coord(inedge(iflowin,1),:))';
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%FUNCTION "getsurnode"
%--------------------------------------------------------------------------
function [esurn,nsurn] = getsurnode(inode,esurn1,esurn2,nsurn1,nsurn2)
%"countelem" counts how many elements there are surrounding each node 
%evaluated.
countelem = esurn2(inode + 1) - esurn2(inode);
%"esurn" calculates what are the elements surrounding each node evaluated.
ie = 1:countelem;
esurn(ie) = esurn1(esurn2(inode) + ie);

%"countnode" counts how many nodes there are surrounding each node
%evaluated.
countnode = nsurn2(inode + 1) - nsurn2(inode);
%"nsurn" calculates what are the nodes surrounding each node 
%evaluated.
in = 1:countnode;
nsurn(in) = nsurn1(nsurn2(inode) + in);

%--------------------------------------------------------------------------
%FUNCTION "shiftchoosen"
%--------------------------------------------------------------------------
function [vecout] = shiftchoosen(vecin,numelem,letter)
%Initialize "vecout". It is the vector reordered according vector element
%choosen.
vecout(1:length(vecin),1) = vecin;

%It finds the position of element choosen in vector "vecin"
if strcmp(letter,'val')
    %Define "i" (auxiliary conuter 1:length(vecin))
    i = 1:length(vecin);
    %Poits its position
    pointpos = i(logical(vecin == numelem));
%"pointpos" receives the position
else
    pointpos = numelem;
end  %End of IF

%It uses "circshift" to reorder the vector
circpos = length(vecin) - pointpos + 1;
vecout = circshift(vecout,circpos);

%--------------------------------------------------------------------------
%FUNCTION "getcentelem"
%--------------------------------------------------------------------------
%This function calculate the coordinate of center's element evaluated 
%(xunkn, yunkn, zunkn). A vector "coordunkn" (3x1) is returned.
%"coordunkn(1)" does mean the coordinate x of unknown point.
%The inflow data is "elemnode" which is a vector with all nodes which
%constitute the element evaluated.
%That is the coordinate of each node who constitute the element evaluated.
%"elemnnode(4)" is the value of fourth column of "elem" matrix. If the elem 
%is a triangle "elemnode(4)" is 0 and the statment (elemnode(4) > 0) is 0. 
%If the elem is a quadrangle the statment is 1 and is added to 3. 
function [centelem] = getcentelem(coord,elem)                         
%Initialize the matrix "centelem"
centelem = zeros(size(elem,1),3);

%Define the coordinate to each element
for ielem = 1:(size(elem,1))
    %"elemnode" receives three or four nodes which constitute the element
    elemnode = elem(ielem,1:sum(elem(ielem,1:4) ~= 0));
    %Calculate the centroid
    centelem(ielem,:) = mean(coord(elemnode,:));
end  %End of FOR (each element)

%--------------------------------------------------------------------------
%Function "calcelemarea"
%--------------------------------------------------------------------------
%This function calculate the area of each element. The element is recognaze
%due parameter "numelem"
function [elemarea] = calcelemarea(coord,elem)
%Initialize the vector
elemarea = zeros(size(elem,1),1);
%Swept all elements of domain
for ielem = 1:size(elem,1)
    %Calculate the element's area to quadrangular element
    if elem(ielem,4) > 0
        %"atriang1" calculates the area of first internal triangle
        atriang1 = 0.5*norm(cross((coord(elem(ielem,2),:) - ...
            coord(elem(ielem,1),:)),(coord(elem(ielem,2),:) - ...
            coord(elem(ielem,3),:))));
        %"atriang2" calculates the area of second internal triangle
        atriang2 = 0.5*norm(cross((coord(elem(ielem,4),:) - ...
            coord(elem(ielem,3),:)),(coord(elem(ielem,4),:) - ...
            coord(elem(ielem,1),:))));
        %The element's area is constituted by sum of two calc. areas above. 
        %This case worth whenever the element is a quadrangle 
        elemarea(ielem) = atriang1 + atriang2;
    %Calculate the element's area to triangular element
    elseif (elem(ielem,4) == 0)&&(elem(ielem,3) > 0)
        %"elemarea" is calculated in the one step
        elemarea(ielem) = ...
            0.5*norm(cross((coord(elem(ielem,2),:) - ...
            coord(elem(ielem,1),:)),(coord(elem(ielem,2),:) - ...
            coord(elem(ielem,3),:))));
    elseif (elem(ielem,3) == 0)
        vfrat=coord(elem(ielem,2),:)-coord(elem(ielem,1),:);
        lfrat=norm(vfrat);
        
        fratopt = fopen('fractureoptions.opt','r'); 
        ntypefrat = textscan(fratopt,'%u',1,'headerlines',1);
        ntypefrat=cell2mat(ntypefrat);
        afrat = textscan(fratopt,'%f',ntypefrat,'headerlines',2);
        fclose(fratopt);
        
        afrat=cell2mat(afrat);
        
        ff=elem(ielem,5)-2000;
        
        elemarea(ielem)=afrat(ff)*lfrat;
  
        %-----------------------------------------------------------------%
    end  %End of IF
end  %End of FOR

%--------------------------------------------------------------------------
%Function "getesure"
%--------------------------------------------------------------------------
function [esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,...
    inedge,esurn1,esurn2,nsurn1,nsurn2)
%Initialize the pointers "esureface2" and "esurefull2"
esureface2 = 0;
esurefull2 = 0;
%Initialize counters
mface = 0;
mfull = 0;
%Sewpt all control volumes
for ielem = 1:size(elem,1)
    %Initialize "neighborelemface" and "neighborelemfull"
    neighborelemface = 0;
    neighborelemfull = 0;
    %Get the amount of edge in each control volume
    amountedge = sum(logical(elem(ielem,1:4) ~= 0));
    %Initialize the auxiliary counter
    m = 1;
    cface = 1;
    cfull = 1;
    boundflagface = 0;
    boundflagfull = 0;
    %Define "elemvertices", "vertex" and "next"
    elemvertices = elem(ielem,1:amountedge);
    
    %Swept the internal half-edges
    for j = 1:amountedge + 1
        %In each face, define the first and the second ("next") vertices.
        vertex = elemvertices(1);
        next = elemvertices(2);
        %Get the "inedge" row for edge evaluated
        pointrow = logical(inedge(:,1) == min(vertex,next) & ...
            inedge(:,2) == max(vertex,next));
        %Verify if "pointrow" belongs to "inedge"
        if any(pointrow)
            vprow = find(pointrow~=0);
            for y=1:size(vprow,1)
                prow=vprow(y);
                %Get the neighbor elements:
                %It stores initialy neighbor of faces but put the neighbor of
                %vertices (if "m" == 2, see below)
                oe=inedge(prow,2 + find(inedge(prow,3:4) ~= ielem));
                if size(oe,2)>1
                    oue=min(oe);
                else
                    oue=oe;
                end
                neighborelemfull(cfull) = oue;
                %It stores only neighbor of face.
                neighborelemface(cface) = oue;
                %when "m" is equal to 2, verify if exists vertex neighboring
                if m == 2
                    %Initialize "neighboraux"
                    neighboraux = 0;
                    %Get the elements surrounding the first vertex of second
                    %edge evaluated.
                    [esurn,] = getsurnode(vertex,esurn1,esurn2,nsurn1,nsurn2);
                    %Reorder the elements position in "esurn"
                    [newesurn] = shiftchoosen(esurn,...
                        neighborelemfull(length(neighborelemfull) - 1),'val');
                    %Evaluate if there is vertex neighbor
                    extraelem = setdiff(newesurn,...
                        [ielem neighborelemfull(length(neighborelemfull) - 1:...
                        length(neighborelemfull))]);
                    %There are extra elements
                    if any(extraelem)
                        neighboraux(1) = ...
                            neighborelemfull(length(neighborelemfull) - 1);
                        neighboraux(2:length(extraelem) + 1) = extraelem;
                        neighboraux(length(extraelem) + 2) = ...
                            neighborelemfull(length(neighborelemfull));
                        %Attribute to "neighborelem" the new neighbor
                        neighborelemfull(cfull - 1:cfull + ...
                            length(neighboraux) - 2) = neighboraux;
                        %Redifine "m" and "c"
                        m = 1;
                        cfull = cfull + length(extraelem);
                    end  %End of IF
                end  %End of IF
                %Update "m" and "c"
                m = m + 1;
                cface = cface + 1;
                cfull = cfull + 1;
            end
        %The edge belongs to "bedge"
        else
            %Decrement "m" in order it does not reach 2.
            m = 1;
            boundflagface = cface;
            boundflagfull = cfull;
        end  %End of IF
        
        %It shifts the vertices position in each control volume:
        %The control volume is a Triangle
        if amountedge == 3
            elemvertices = circshift(elemvertices,1:3);
        %The control volume is a Quadrangle
        else
            elemvertices = circshift(elemvertices,2:4);
        end  %End of IF
    end  %End of internal FOR (each control volume)
    %Ensure there is no repeated element
    neighborelemface = unique(neighborelemface,'stable');
    neighborelemfull = unique(neighborelemfull,'stable');

    %Reorder adequadely "neighborelemface"
    if boundflagface ~= 0 && boundflagface <= length(neighborelemface)
        neighborelemface = shiftchoosen(neighborelemface,...
            neighborelemface(boundflagface),'val');
    end  %End of IF ("boundflagface")
    %Reorder adequadely "neighborelemfull"
    if boundflagfull ~= 0 && boundflagfull <= length(neighborelemfull)
        neighborelemfull = shiftchoosen(neighborelemfull,...
            neighborelemfull(boundflagfull),'val');
    end  %End of IF ("boundflagface")
        
    %Attribute the contribution from "esureface" and "esurefull" to vectors 
    %"esureface1" and "esurefull1"
    esureface1(mface + 1:mface + length(neighborelemface)) = ...
        neighborelemface;
    esurefull1(mfull + 1:mfull + length(neighborelemfull)) = ...
        neighborelemfull;
    
    %Attribute the contribution from "esureface" and "esurefull" to vectors 
    %"esureface2" and "esurefull2"
    esureface2 = vertcat(esureface2,length(esureface1));
    esurefull2 = vertcat(esurefull2,length(esurefull1));

    %Increment counters "mface" and "mfull"
    mface = length(esureface1);
    mfull = length(esurefull1);
end  %End of FOR (Swept all elements)

%Turn "esureface1" and "esurefull1" vertical vector
esureface1 = esureface1';
esurefull1 = esurefull1';

%--------------------------------------------------------------------------
%FUNCTION "getrowposition"
%--------------------------------------------------------------------------

function [rowposit] = getrowposition(bedge,inedge,nsurn1,nsurn2)
%Initialize "rowposit"
rowposit = zeros(length(nsurn1),1);
%Initialize "bedgesize"
bedgesize = size(bedge,1);

%Swept the boundary edges
for i = 1:bedgesize
    %Get the vertices:
    vertices = bedge(i,1:2);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(1)) + 1:nsurn2(vertices(1) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    pointvtx = logical(nsurn == vertices(2));
    nsurnpos = nsurnsize(pointvtx);
    %Fill "rowposit"
    rowposit(nsurn2(vertices(1)) + nsurnpos) = i;

    %----------------------
    %Evaluate the vertex 2:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(2)) + 1:nsurn2(vertices(2) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    pointvtx = logical(nsurn == vertices(1));
    nsurnpos = nsurnsize(pointvtx);
    %Fill "rowposit"
    rowposit(nsurn2(vertices(2)) + nsurnpos) = i;
end  %End of FOR ("bedge")

%Swept the boundary edges
for i = 1:size(inedge,1)
    %Get the vertices:
    vertices = inedge(i,1:2);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(1)) + 1:nsurn2(vertices(1) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    pointvtx = logical(nsurn == vertices(2));
    nsurnpos = nsurnsize(pointvtx);
    %Fill "rowposit"
    rowposit(nsurn2(vertices(1)) + nsurnpos) = i + bedgesize;

    %----------------------
    %Evaluate the vertex 2:
    
    %Get the nsurn for the vertex evaluated.
    nsurn = nsurn1(nsurn2(vertices(2)) + 1:nsurn2(vertices(2) + 1));
    %Get a vector with the size of "nsurn"
    nsurnsize = 1:length(nsurn);
    %Identify the vertex 2 position in "nsurn"
    pointvtx = logical(nsurn == vertices(1));
    nsurnpos = nsurnsize(pointvtx);
    %Fill "rowposit"
    rowposit(nsurn2(vertices(2)) + nsurnpos) = i + bedgesize;
end  %End of FOR ("bedge")


%-------------------------------------------------------------------------%
%FUNCTION "reordernsurn"
%-------------------------------------------------------------------------%
function [nsurn1,nsurn2] = reordernsurn(esurn1,esurn2,...
            nsurn1,nsurn2,bedge,coord,elem)
y=1;
fratnodes=zeros(1,1);
numfrat=sum(elem(:,5)>2000);
for i=1:numfrat
    for j=1:2
        a=elem(i,j);
        if isempty(find(fratnodes==a))==1
            fratnodes(y)=a;
            y=y+1;
        end
    end
end

for i=1:size(fratnodes,2)
    a=fratnodes(i);
    if a~=0
        va=nsurn1(nsurn2(a)+1:nsurn2(a+1));
        c=0;
        for y=1:size(va,2)
            if isempty(find(fratnodes==va(y)))==0
                vc(c+f)=va(y);
                c=c+1;
            end
        end
        for g=size(nsurn1,2):-1:nsurn2(a)+1
            nsurn1(g+c)=nsurn1(g);
        end
        for g=a+1:size(nsurn2,1)
            nsurn2(g)=nsurn2(g)+c;
        end
        for s=1:c
            nsurn1(nsurn2(a)+1+s)=vc(s);
        end
    end
end

nsurn1 = nsurn1'; 

is_bound = zeros(size(coord,1),1);

for i=1:size(bedge,1),
   is_bound(bedge(i,1))=1; 
end

for i=1:size(coord,1),
    if is_bound(i)==0
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if (elem(prim,4)==0)&&(elem(prim,3)==0)
            bp=2;
        elseif (elem(prim,4)==0)&&(elem(prim,3)~=0)
            bp=3;
        else
            bp=4;
        end
        if (elem(ult,4)==0)&&(elem(ult,3)==0)
            bu=2;
        elseif (elem(ult,4)==0)&&(elem(ult,3)~=0)
            bu=3;
        else
            bu=4;
        end
        for j=1:bp,
            for u=1:bu,
                if ((elem(ult,u)==elem(prim,j))&&(elem(prim,j)~=i))
                    nnsn(1)=elem(prim,j);
                end
            end
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if (elem(atual,4)==0)&&(elem(atual,3)==0)
                bp=2;
            elseif (elem(atual,4)==0)&&(elem(atual,3)~=0)
                bp=3;
            else
                bp=4;
            end
            if (elem(ant,4)==0)&&(elem(ant,3)==0)
                bu=2;
            elseif (elem(ant,4)==0)&&(elem(ant,3)~=0)
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end     
    else 
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if (elem(prim,4)==0)&&(elem(prim,3)==0)
            bp=2;
        elseif (elem(prim,4)==0)&&(elem(prim,3)~=0)
            bp=3;
        else
            bp=4;
        end
        if (elem(ult,4)==0)&&(elem(ult,3)==0)
            bu=2;
        elseif (elem(ult,4)==0)&&(elem(ult,3)~=0)
            bu=3;
        else
            bu=4;
        end
        for j=1:size(bedge,1)
            for p=1:bp,
                if (i==elem(prim,p))
                    if p==bp
                        n1=elem(prim,1);
                    else
                        n1=elem(prim,p+1);
                    end
                end
            end
            nnsn(1)=n1;            
            for u=1:bu,
                if (i==elem(ult,u))
                    if u==1
                        n1=elem(ult,bu);
                    else
                        n1=elem(ult,u-1);
                    end
                end
            end
            nnsn((nsurn2(i+1)-nsurn2(i)))=n1;
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if (elem(atual,4)==0)&&(elem(atual,3)==0)
                bp=2;
            elseif (elem(atual,4)==0)&&(elem(atual,3)~=0)
                bp=3;
            else
                bp=4;
            end
            if (elem(ant,4)==0)&&(elem(ant,3)==0)
                bu=2;
            elseif (elem(ant,4)==0)&&(elem(ant,3)~=0)
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end
    end            
end
