function regProduct = normcorre_registration(wholeSession)
%========================================================================================================================================
%
% Uses the NoRMCorre algorithm to apply a rigid motion correction to pre-processed imaging data. It must be saved as a .mat file 
% containing the variables 'trialType', 'origFileNames', 'expDate', and either 'wholeSession' or 'regProduct' (the latter in case I want
% to register files that have already been processed by my old registration algorithm). The imaging data array (which must have 
% dimensions [y, x, plane, volume, trial])is the only variable that is actually used, the rest are just re-saved as a new .mat file in 
% the same directory as the source file.
% 
%=========================================================================================================================================

% Set parameters % Determined by MM
options_rigid = NoRMCorreSetParms('d1', size(wholeSession, 1), 'd2', size(wholeSession, 2), 'd3', size(wholeSession, 3), ...
                    'max_shift', [25, 25, 2], ...
                    'init_batch', 100, ...
                    'us_fac', 50 ...
                    ); 
% options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [20, 20, 2], ...
%                     'init_batch', 100, ...
%                     'grid_size', [100, 100] ...
%                     );
                
% Rigid registration
tic; [regProduct, ~, ~, ~] = normcorre(wholeSession, options_rigid); toc

end