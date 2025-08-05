function func_Preproc(inpfiles, Dirs, subj, ses, AcqParams)

save_path = Dirs.out; 
subj_name = subj.id;
TR = AcqParams.tr; 

%% Run Preprocessing Batch
spm fmri;

matlabbatch = func_PreprocBatch(inpfiles, AcqParams, Dirs);
spm_jobman('run', matlabbatch)

%% Deleting unnecessary files and moving results to the related folder
 
% Functional data
funcSaveDir = fullfile(save_path, subj_name, ses, ...
    'func');
mkdir(funcSaveDir);
datapath = fullfile(subj.path, ses, 'func');
delete(fullfile(datapath, 'a*.*'));
delete(fullfile(datapath, 'meana*.*'));
delete(fullfile(datapath, 'smeana*.*'));
delete(fullfile(datapath, '*.mat'));
movefile(fullfile(datapath, 'sra*.*'), funcSaveDir);
movefile(fullfile(datapath, 'rp_*.txt'), funcSaveDir);
movefile(fullfile(datapath, 'ra*.nii'), funcSaveDir);

 
%% Run Artifact Detection Batch
  
clear matlabbatch;
[matlabbatch, art_pth] = func_ArtDetection_batch(funcSaveDir, save_path, subj_name, ses, TR);
spm_jobman('serial', matlabbatch);
art_batch(fullfile(art_pth, 'SPM.mat'));