function CopiaDeSeguranca( S_old,VPI,countime,step,watercut,oilrecovery,cumulateoil,time )

global problemData

fid = fopen([char(problemData.problem.outputPath),filesep,'CopiaDeSeguranca.txt'],'w');

%Write data related to VPI
fprintf(fid,'VPI \r\n');
fprintf(fid,'%26.16E \r\n',VPI);
%Jump a line

%Write data related to Countime
fprintf(fid,'COUNTIME \r\n');
fprintf(fid,'%26.16E \r\n',countime);
%Jump a line

%Write data related to Step
fprintf(fid,'STEP \r\n');
fprintf(fid,'%26.16E \r\n',step);
%Jump a line

%Write data related to Saturation
fprintf(fid,'SATURATION \r\n');
fprintf(fid,'%26.16E \r\n',S_old);
%Jump a line

%Write data related to watercut
fprintf(fid,'WATERCUT \r\n');
fprintf(fid,'%26.16E \r\n',watercut);
%Jump a line

%Write data related to oilrecovery
fprintf(fid,'OILRECOVERY \r\n');
fprintf(fid,'%26.16E \r\n',oilrecovery);
%Jump a line

%Write data related to cumulateoil
fprintf(fid,'CUMULATEOIL \r\n');
fprintf(fid,'%26.16E \r\n',cumulateoil);
%Jump a line

%Write data related to time
fprintf(fid,'TIME \r\n');
fprintf(fid,'%26.16E \r\n',time);
%Jump a line

fclose(fid);

end

