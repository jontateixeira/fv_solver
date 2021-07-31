function [ S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time ] = ...
         condinicial( S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time )
%
global problemData
res_folder = char(problemData.problem.outputPath);
results = [res_folder,filesep,'Saturation.mat']; 
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'VPI.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'Countime.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'Oilrecovery.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'Cumulateoil.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'Watercut.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'TimeStep.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'DT.mat'];
if exist(results,'file') ~= 0
    load(results);
end

results = [res_folder,filesep,'TIME.mat'];
if exist(results,'file') ~= 0
    load(results);
end

end

