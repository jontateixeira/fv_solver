%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 01/05/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza; Luiz Eduardo
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION  
%--------------------------------------------------------------------------

function postprocessor(passo,timestep,coord,elem,pressure,watersaturation)
%
global problemData eps_

%--------------------------------------------------------------------------
%Create the post file (*.vtk) to PRESSURE Field 
%--------------------------------------------------------------------------
%Initialization parameters

%Attribution of file extension
ext = '.vtk';
%First file name
fname = 'PressOut';

%--------------------------------------------------------------------------
%Structure's parameters

%Number of nodes to used
nnode = size(coord,1);
%"countriang" count how many triangle there is in the domain
countriang = sum(elem(:,4) == 0);
%"countquad" count how many quadrangles there is in the domain
countquad = sum(elem(:,4) ~= 0);
%Initialize "type_elem"
type_elem = zeros(countriang+countquad,1);
%Number of elements to be used
nelem = size(type_elem,1);
%"sumnodebyelem" is a parameter which verify how many nodes constitute each
%element multiplied by element amount
sumnodebyelem = sum(sum(elem(1:nelem,1:4) ~= 0));
pressure=pressure(1:nelem)*1e4; % Conversão para 'bar'!!!!!!!
watersaturation=watersaturation(1:nelem);

%--------------------------------------------------------------------------
%Write the file

%The command below create the file name joining the parameters which follow:
%"fname" (file name), a tag ("00"), "step" (used in case where is wanted 
%one result by timestep) and the file extension ("step") 
fname_vtk = [fname, sprintf('%05d',passo), ext];
fprintf('** creating %s file\n', fname_vtk);
%Open the file to write the text and create the writer "fid".
%The use of letter "w" implies in a quite write, erasing any word in a
%existent file
fid = fopen([char(problemData.problem.outputPath),filesep,fname_vtk],'w'); 
%Write head informations
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data \r\n');
fprintf(fid,'ASCII \r\n');
%Information about grid type
fprintf(fid,'DATASET UNSTRUCTURED_GRID \r\n');
fprintf(fid,'FIELD FieldData 2\n');
fprintf(fid,'TIME 1 1 float\n');
fprintf(fid,'%f\n',timestep);
fprintf(fid,'CYCLE 1 1 int\n');
fprintf(fid,'%d\n',passo);

%Write the POINT informations:
%Head (POINT)
fprintf(fid,'POINTS %i float \r\n\r\n',nnode);
%Distribution (POINT)
data1 = [coord(:,1:3)'];
%Print the distribution
fprintf(fid,'%26.16E %26.16E %26.16E \r\n',data1);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL informations:
%Head (CELL)
fprintf(fid,'CELLS %i %i \r\n\r\n',nelem,sumnodebyelem + nelem);

%There is trianglular elements
if countriang > 0
    %Define how many nodes constitute each element
    %"nodebyelem" is a parameter which verify how many nodes constitute 
    %each element
    nodebyelem = 3;
    %"numdata" is the first column of data to be printed in the CELL sect.
    numdata = nodebyelem*ones(countriang,1);

    %Considering which the initial data (counter) begining from zero, the 
    %node number is diminish of one (-1) in the node definition (command 
    %below). The decision command below is used in order to define if the 
    %elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
    %("nodebyelem" == 4)
    data2 = ...
        [numdata';(elem(1:countriang,1) - 1)';...
        (elem(1:countriang,2) - 1)';...
        (elem(1:countriang,3) - 1)'];
        fprintf(fid,'%i %i %i %i \r\n',data2);
    %Definition of element type (5 is a code to triangle)
    type_elem(1:countriang) = 5;
end  %End of IF (triangles)

%There is quadrangle element
if countquad > 0
    %Define how many nodes constitute each element
    %"nodebyelem" is a parameter which verify how many nodes constitute 
    %each element
    nodebyelem = 4;
    %"numdata" is the first column of data to be printed in the CELL sect.
    numdata = nodebyelem*ones(countquad,1);

    %Considering which the initial data (counter) begining from zero, the 
    %node number is diminish of one (-1) in the node definition (command 
    %below). The decision command below is used in order to define if the 
    %elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
    %("nodebyelem" == 4)
    data2 = ...
        [numdata';(elem(countriang + 1:countriang+countquad,1) - 1)';...
        (elem(countriang+1:countriang+countquad,2) - 1)';...
        (elem(countriang+1:countriang+countquad,3) - 1)';...
        (elem(countriang+1:countriang+countquad,4) - 1)'];
    fprintf(fid,'%i %i %i %i %i \r\n',data2);
    %Definition of element type (7 is a code to quadrangle)
    type_elem(countriang + 1:nelem) = 7;            
end  %End of IF (quadrangle)

p=pressure(1:nelem,1);
S=watersaturation(1:nelem); S(S < eps_) = 0;

%Jump a line
fprintf(fid,'\r\n');

%Cell type information
fprintf(fid,'CELL_TYPES %i \r\n\r\n',nelem);
fprintf(fid,'%i \r\n',type_elem);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL_DATA informations:
%Head (CELL_DATA)
fprintf(fid,'CELL_DATA %i \r\n',nelem);
%Write data related to PRESSURE
fprintf(fid,'SCALARS PRESSURE float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',full(p));
%Jump a line
fprintf(fid,'\r\n');

%fprintf(fid,'CELL_DATA %i \r\n',nelem);
%Write data related to SATURATION
fprintf(fid,'SCALARS SATURATION float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',S);
%Jump a line

fclose(fid);

end  %End of function

