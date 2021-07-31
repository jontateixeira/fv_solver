function startup_SOLVER
%Amend MATLAB PATH to handle functions implementation.

% $Date: 2020-07-02$
% $Revision: 1 $

   % current folder
   ROOT_DIR = fileparts(mfilename('fullpath'));
   
   % add common folder and subfolders(genpath)
   addpath(genpath(fullfile(ROOT_DIR,'preprocess')));
   
   % add framework folder and subfolders(genpath)
   addpath(genpath(fullfile(ROOT_DIR,'framework')));
   
   % add postprocess folder and subfolders(genpath)
   addpath(genpath(fullfile(ROOT_DIR,'postprocess')));
   
   % add util folder and subfolders(genpath)
   addpath(genpath(fullfile(ROOT_DIR,'utils')));
   
   % display helloworld mensage
   disp('* FV_SOLVER were set into MATLAB PATH')
end
