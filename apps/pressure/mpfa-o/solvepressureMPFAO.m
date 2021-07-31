function [pressure,influx,bflux,q] = solvepressureMPFAO(transmvecleft,...
      transmvecright,knownvecleft,knownvecright,storeinv,Bleft,Bright,...
      wells,mapinv,maptransm,mapknownvec,pointedge,mobility,bodyterm)
%Define global parameters:
global elem bedge inedge phflw

if phflw==0 || phflw==1
   phasekey=1;
elseif phflw==2
   phasekey=2; 
end

%Initialize the global matrix which have order equal to parameter 
%"size(elem)".
M = sparse(size(elem,1),size(elem,1));
%Initialize "mvector" which is the independent vector of algebric system.
mvector = sparse(size(elem,1),1);

%Initialize "updtmatcoeffleft". It stores the matrices coefficients in 
%order to calculate the velocity.
updtmatcoeffleft = zeros(length(transmvecleft),1);
%"updtvectorleft" receives the known terms
updtvectorleft = zeros(length(knownvecleft),1);
%Initialize counters
countcoeff = 0;
countvec = 1;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%--------------------------------------------------------------------------
%Swept the edges on the boundary and domain

%Swept "bedge"
for i = 1:bedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = bedge(i,1:2);
    %Get the element on the left
    leftelem = bedge(i,3);
    
    %----------------------------------------------------------------------
    %Vertex 1:
    
    %Get the transmissibility coefficients:
    [mtxleft,vecleft,esurn] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,i,1,mobility,bedgesize,...
        inedgesize);

    %-----------------------------
    %Assembly the algebraic system

    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;

    %-----------------------------
    
    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;

    %----------------------------------------------------------------------
    %Vertex 2:
    
    %Get the transmissibility coefficients:
    [mtxleft,vecleft,esurn] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,i,2,mobility,bedgesize,...
        inedgesize);

    %-----------------------------
    %Assembly the algebraic system

    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;
end  %End of FOR (Swept boundary edges)

%Swept "inedge"
for i = 1:inedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = inedge(i,1:2);
    %Get the element on the left
    leftelem = inedge(i,3);
    %Get the element on the right
    rightelem = inedge(i,4);
    
    %----------------------------------------------------------------------
    %Vertex 1:
    
    %Get the transmissibility coefficients (left contribution):
    [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,bedgesize + i,1,...
        mobility,bedgesize,inedgesize);

    %Get the transmissibility coefficients (right contribution):
    [mtxright,vecright,esurn] = gettransmcoeff(transmvecright,...
        knownvecright,storeinv,Bright,2,mapinv,maptransm,mapknownvec,...
        pointedge,bodyterm,vertices(1),vertices(2),phasekey,...
        bedgesize + i,1,mobility,bedgesize,inedgesize);

    %-----------------------------
    %Assembly the algebraic system
    
    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    M(rightelem,esurn) = M(rightelem,esurn) + mtxright';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;
    mvector(rightelem) = mvector(rightelem) - vecright;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;

    %----------------------------------------------------------------------
    %Vertex 2:
    
    %Get the transmissibility coefficients (left contribution):
    [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,bedgesize + i,2,...
        mobility,bedgesize,inedgesize);

    %Get the transmissibility coefficients (right contribution):
    [mtxright,vecright,esurn] = gettransmcoeff(transmvecright,...
        knownvecright,storeinv,Bright,2,mapinv,maptransm,mapknownvec,...
        pointedge,bodyterm,vertices(2),vertices(1),phasekey,...
        bedgesize + i,2,mobility,bedgesize,inedgesize);
    
    %-----------------------------
    %Assembly the algebraic system
    
    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    M(rightelem,esurn) = M(rightelem,esurn) + mtxright';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;
    mvector(rightelem) = mvector(rightelem) - vecright;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;
end  %End of FOR (Swept internal edges)

%Release Memory
clear transmvecleft transmvecright knownvecleft knownvecright storeinv ...
    Bleft Bright;

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,mvector] = addsource(sparse(M),mvector,wells);

%% problema numcase= 15.7; % gao e Wu 2013 e Terekhov 2016
%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver". 
%This function returns the pressure field with value put in each colocation 
%point.
[pressure] = solver(M,mvector);

%Calculate flow rate through edge. "satkey" equal to "1" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[influx,bflux,q] = calcflowrateMPFA(updtmatcoeffleft,updtvectorleft,pressure);


