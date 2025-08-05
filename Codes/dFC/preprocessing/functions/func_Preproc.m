function func_Preproc(inpfiles, Dirs, subj, ses, AcqParams)

save_path = Dirs.out; 
subj_name = subj.name;
TR = AcqParams.tr; 

%% Run Preprocessing Batch
spm fmri;

matlabbatch = func_PreprocBatch(inpfiles, AcqParams, Dirs);
spm_jobman('run', matlabbatch)

%% Deleting unnecessary files and moving results to the related folder
 
% Functional data
motionCorrectedDir = fullfile(save_path, subj_name, ses, ...
    'func');
mkdir(motionCorrectedDir);
datapath = fullfile(subj.dir, ses, 'func');
delete(fullfile(datapath, 'a*.*'));
delete(fullfile(datapath, 'meana*.*'));
delete(fullfile(datapath, 'wmeana*.*'));
delete(fullfile(datapath, 'swmeana*.*'));
delete(fullfile(datapath, '*.mat'));
delete(fullfile(datapath, 'ra*.*'));
movefile(fullfile(datapath, 'swra*.*'), motionCorrectedDir);
movefile(fullfile(datapath, 'rp_*.txt'), motionCorrectedDir);
movefile(fullfile(datapath, 'wra*.*'), motionCorrectedDir);
 
% Structural Data
stresDir = fullfile(save_path, subj_name, ses, 'anat');
mkdir(stresDir);
datapath = fullfile(subj.dir, ses, 'anat');
movefile(fullfile(datapath, 'mri'), stresDir);
movefile(fullfile(datapath, 'report'), stresDir);
 
%% Run Artifact Detection Batch
  
clear matlabbatch;
[matlabbatch, art_pth] = func_ArtDetection_batch(motionCorrectedDir, save_path, subj_name, ses, TR);
spm_jobman('serial', matlabbatch);
art_batch(fullfile(art_pth, 'SPM.mat'));