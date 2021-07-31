%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: %This funtion returns parameters surrounding a node evaluated 
%(elements and nodes)
  

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [esurn,nsurn,fsurn] = getsurnode(inode)
%Define global parameters:
global esurn1 esurn2 nsurn1 nsurn2 keyfrac;

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

if strcmp(keyfrac,'y')%BSB-i
    %Add global parameters:
    global fsurn1;
    
    %"countnode" counts how many nodes there are surrounding each node
    %evaluated.
    countnode = nsurn2(inode + 1) - nsurn2(inode);
    %"fsurn" calculates what are the fractures surrounding each node
    %evaluated.
    in = 1:countnode;
    fsurn(in) =fsurn1(nsurn2(inode) + in);
else
    fsurn=zeros(length(nsurn),1);
end
