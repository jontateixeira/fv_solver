function [influx,bflux,q] = calcflowrateMPFA(updatecoeffleft,updatevectorleft,pressure)
%Define global parameters:
global elem bedge inedge

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "counters"
countmtrx = 0;
countvec = 1;

%Initialise "flowresult". This accumulates the flowrate in all elements 
%associated to producer wells. It's a flow rate resulting for each element. 
q = zeros(size(elem,1),1);

%--------------------------------------------------------------------------
%Switch the way to calculate "flowrate" according to "phflw" and
%"satkey" values

%--------------------------------------------------------------------------
influx = zeros(inedgesize,1);
bflux = zeros(bedgesize,1);

%Swept the boundary edges
for i = 1:bedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = bedge(i,1:2);
    %Define the element on the left
    elemleft = bedge(i,3);

    %------------------------------------------------------------------
    %Vertex 1:

    %Get the elements and the nodes surrounding the vertex 1
    [esurn,] = getsurnode(vertices(1));

    %Get "mtxleft" and "vecleft"
    mtxleft = updatecoeffleft(countmtrx + 1:countmtrx + length(esurn));
    vecleft = updatevectorleft(countvec);

    %Calculate the flowrate for edge evaluated
    getflowrate = mtxleft'*pressure(esurn) + vecleft;
    %Attribute to "flowrate" the darcy velocity
    bflux(i) = bflux(i) + getflowrate;  

    %Attribute to "flowresult" the darcy velocity
    q(elemleft) = q(elemleft) + getflowrate;  

    %Update the counters:
    countmtrx = countmtrx + length(esurn);
    countvec = countvec + 1;

    %------------------------------------------------------------------
    %Vertex 2:

    %Get the elements and the nodes surrounding the vertex 1
    [esurn,] = getsurnode(vertices(2));

    %Get "mtxleft" and "vecleft"
    mtxleft = updatecoeffleft(countmtrx + 1:countmtrx + length(esurn));
    vecleft = updatevectorleft(countvec);

    %Calculate the flowrate for edge evaluated
    getflowrate = mtxleft'*pressure(esurn) + vecleft;
    %Attribute to "flowrate" the darcy velocity
    bflux(i) = bflux(i) + getflowrate;  

    %Attribute to "flowresult" the darcy velocity
    q(elemleft) = q(elemleft) + getflowrate;  

    %Update the counters:
    countmtrx = countmtrx + length(esurn);
    countvec = countvec + 1;
end  %End of FOR (Swept boundary edges)

%Swept the internal edges
for i = 1:inedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = inedge(i,1:2);
    %Define the elements on the left and on the right
    elemleft = inedge(i,3);
    elemright = inedge(i,4);

    %------------------------------------------------------------------
    %Vertex 1:

    %Get the elements and the nodes surrounding the vertex 1
    [esurn,] = getsurnode(vertices(1));

    %Get "mtxleft" and "vecleft"
    mtxleft = updatecoeffleft(countmtrx + 1:countmtrx + length(esurn));
    vecleft = updatevectorleft(countvec);

    %Calculate the flowrate for edge evaluated
    getflowrate = mtxleft'*pressure(esurn) + vecleft;
    %Attribute to "flowrate" the darcy velocity
    influx(i) = influx(i) + getflowrate;  

    %Attribute to "flowresult" the darcy velocity (on the left)
    q(elemleft) = q(elemleft) + getflowrate;  
    %Attribute to "flowresult" the darcy velocity (on the right)
    q(elemright) = q(elemright) - getflowrate;  

    %Update the counters:
    countmtrx = countmtrx + length(esurn);
    countvec = countvec + 1;

    %------------------------------------------------------------------
    %Vertex 2:

    %Get the elements and the nodes surrounding the vertex 1
    [esurn,] = getsurnode(vertices(2));

    %Get "mtxleft" and "vecleft"
    mtxleft = updatecoeffleft(countmtrx + 1:countmtrx + length(esurn));
    vecleft = updatevectorleft(countvec);

    %Calculate the flowrate for edge evaluated
    getflowrate = mtxleft'*pressure(esurn) + vecleft;
    %Attribute to "flowrate" the darcy velocity
    influx(i) = influx(i) + getflowrate;  

    %Attribute to "flowresult" the darcy velocity (on the left)
    q(elemleft) = q(elemleft) + getflowrate;  
    %Attribute to "flowresult" the darcy velocity (on the right)
    q(elemright) = q(elemright) - getflowrate;  

    %Update the counters:
    countmtrx = countmtrx + length(esurn);
    countvec = countvec + 1;
end  %End of FOR ("inedge")
