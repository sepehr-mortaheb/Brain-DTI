clc 
clear 

%% Initialization 

% --- Set the following directories --- 

% Directory of the BIDS formated data:
bids_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/data/raw_bids';
% Save directory of the fMRI processing:
save_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/data/preprocessed';

%% 
% --- Set the Participants and Sessions Information --- 

subjectDirs = dir(fullfile(bids_dir, 'sub-*'));
Subjects = struct();
for i=1:length(subjectDirs)
    SubjName = subjectDirs(i).name;
    SubjPath = fullfile(bids_dir, SubjName);
    
    % Get session folders
    sessionDirs = dir(fullfile(SubjPath, 'ses-*'));
    validSessions = {};
    for j=1:length(sessionDirs)
        sessName = sessionDirs(j).name;
        validSessions{end+1} = sessName;
    end
    Subjects(i).id = SubjName;
    Subjects(i).path = SubjPath;
    Subjects(i).sessions = validSessions; 
end

%% 

for i=1:length(Subjects)
    subj = Subjects(i).id; 
    subj_path = Subjects(i).path;
    n_sessions = numel(Subjects(i).sessions); 
    
    anat_files = cell(n_sessions, 1);
    for s = 1: n_sessions
        ses = Subjects(i).sessions{s};
        anat_files{s} = fullfile(subj_path, ses, 'anat', strcat(subj, '_', ses, '_acq-ge_run-fmri_T1w.nii'));
        if ~exist(anat_files{s}, 'file')
            anat_files{s} = fullfile(subj_path, ses, 'anat', strcat(subj, '_', ses, '_acq-ge_run-dti_T1w.nii'));
        end
    end
    
    if n_sessions == 2
        times = [0 0.5];
    else 
        times = [0 0.5 1];
    end
    
            
    matlabbatch = {};

    matlabbatch{1}.spm.tools.longit.series.vols = anat_files;
    matlabbatch{1}.spm.tools.longit.series.times = times;
    matlabbatch{1}.spm.tools.longit.series.noise = NaN;
    matlabbatch{1}.spm.tools.longit.series.wparam = [0 0 100 25 100];
    matlabbatch{1}.spm.tools.longit.series.bparam = 1000000;
    matlabbatch{1}.spm.tools.longit.series.write_avg = 1;
    matlabbatch{1}.spm.tools.longit.series.write_jac = 0;
    matlabbatch{1}.spm.tools.longit.series.write_div = 0;
    matlabbatch{1}.spm.tools.longit.series.write_def = 0;
    
    spm_jobman('run', matlabbatch)
    
    SaveDir = fullfile(save_dir, subj, 'anat_template');
    if ~exist(SaveDir, 'dir')
        mkdir(SaveDir);
    end
    datapath = fullfile(subj_path, Subjects(i).sessions{1}, 'anat');
    movefile(fullfile(datapath, 'avg_*.*'), SaveDir);
    
    clear matlabbatch;   
end