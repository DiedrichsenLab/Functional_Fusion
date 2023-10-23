function spmj_realign(data,varargin)
% spmj_realign(data)
% INPUT: 
%   data:    cell array (runs) of cell array (images) 
% VARARGINOPTIONS: 
% 'weight',{'weighting image'}

weight=''; 
vararginoptions(varargin,{'weight'}); 

%_______DEFAULTS_________________________________
J.data=data;
J.eoptions.quality = 0.9;                                                                            
J.eoptions.sep = 4;                                                                                  
J.eoptions.fwhm = 5;                                                                                 
J.eoptions.rtm = 0;  % 0, register to the first / 1, register to the mean                                                                                
J.eoptions.interp = 2;                                                                               
J.eoptions.wrap = [0 0 0];                                                                           
J.eoptions.weight = {weight};                                                                            
J.roptions.which = [2 1];                                                                            
J.roptions.interp = 4;                                                                               
J.roptions.wrap = [0 0 0];                                                                           
J.roptions.mask = 1;                                                                                 
J.roptions.prefix = 'r';                                                                             

matlabbatch{1}.spm.spatial.realign.estwrite= J;
spm_jobman('run',matlabbatch);