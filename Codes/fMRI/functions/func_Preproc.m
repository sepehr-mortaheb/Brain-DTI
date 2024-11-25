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
    'func', 'MotionCorrected');
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
delete(fullfile(datapath, 'cat*.mat'));
delete(fullfile(datapath, 'cat*.jpg'));
delete(fullfile(datapath, 'cat*.txt'));
delete(fullfile(datapath, 'cat*.xml'));
delete(fullfile(datapath, 'msub*.*'));
delete(fullfile(datapath, 'wp1*.nii'));
delete(fullfile(datapath, 'wp2*.nii'));
movefile(fullfile(datapath, 'mwp1*.*'), stresDir);
movefile(fullfile(datapath, 'mwp2*.*'), stresDir);
movefile(fullfile(datapath, 'mwp3*.*'), stresDir);
movefile(fullfile(datapath, 'wmsub*.*'), stresDir);
movefile(fullfile(datapath, 'y_*.*'), stresDir);
movefile(fullfile(datapath, 'iy_*.*'), stresDir);
movefile(fullfile(datapath, 'cat*.pdf'), stresDir);
 
%% Run Artifact Detection Batch
  
clear matlabbatch;
[matlabbatch, art_pth] = func_ArtDetection_batch(motionCorrectedDir, save_path, subj_name, ses, TR);
spm_jobman('serial', matlabbatch);
art_batch(fullfile(art_pth, 'SPM.mat'));