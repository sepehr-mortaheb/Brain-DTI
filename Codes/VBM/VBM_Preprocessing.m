%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script for Voxel-Based Morphometry (VBM) analysis for         
% LONGITUDINAL data using CAT12. It contains Segmentation, Normalization, 
% Smoothing. It considers that the data are saved in BIDS format. 
% 
% NOTE: Before running the script, make sure that the cat12 package is open
% and it is set to the Expert Mode.
% 
% Written by: Sepehr Mortaheb 
% Lab for Equilibrium Investigations and Aerospace (LEIA)
% University of Antwerp
% April 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close 

%% SET Directories
spm_addr = '/Users/sjillings/spm12';
data_dir = '/Users/sjillings/Desktop/sepehr/projects/Brain_DTI_VBM/data';
res_dir = '/Users/sjillings/Desktop/sepehr/projects/Brain_DTI_VBM/results2'; 
tmplt_dir = '/Users/sjillings/Desktop/sepehr/projects/Brain_DTI_VBM/Templates';

seg_dir = fullfile(res_dir, 'segmentation'); 
if ~isfolder(seg_dir)
    mkdir(seg_dir)
end

%% Segmentation 

% Reading the Data (for now just modes 1 qnd 3 are supported)
mode = 3; % Mode of data>> 1:one T1 per session for each subject
          %                2:all the T1 images that exist
          %                3:a specific participant 
data = VBM_read_data(mode, data_dir, seg_dir);

% Creating the Segmentation batch 
matlabbatch = VBM_Segmentation_batch(data, spm_addr, tmplt_dir);

% Run Segmentation Batch
spm_jobman('run', matlabbatch)


%% Smoothing 

% Reading the data 
subj_list = dir(seg_dir);
subj_list = subj_list(contains({subj_list.name}, 'sub'));
gm_data = cell(1, 1);
wm_data = cell(1, 1);
csf_data = cell(1, 1);

c = 1;
for i = 1:numel(subj_list)
    subj_dir = fullfile(subj_list(i).folder, subj_list(i).name);
    sessions = dir(subj_dir);
    sessions = sessions(contains({sessions.name}, 'ses'));
    for j = 1:numel(sessions)
        session_dir = fullfile(sessions(j).folder, sessions(j).name);
        anat_dir = fullfile(session_dir, 'anat');
        all_results= dir(anat_dir);
        gm_str = all_results(startsWith({all_results.name}, 'mwmwp1'));
        gm = fullfile(gm_str.folder, gm_str.name);
        gm_data{c, 1} = gm;
        
        wm_str = all_results(startsWith({all_results.name}, 'mwmwp2'));
        wm = fullfile(wm_str.folder, wm_str.name);
        wm_data{c, 1} = wm;
        
        csf_str = all_results(startsWith({all_results.name}, 'mwmwp3'));
        csf = fullfile(csf_str.folder, csf_str.name);
        csf_data{c, 1} = csf;
        c = c+1;
    end
end

% Creating the smoothing batch 
matlabbatch = VBM_Smoothing_batch(vertcat(gm_data, wm_data, csf_data));

% Running the Smoothing batch
spm_jobman('run', matlabbatch)

% Moving the results to the smoothing results folder 
smth_dir = fullfile(res_dir, 'smoothing'); 
if ~isfolder(smth_dir)
    mkdir(smth_dir)
end


for i = 1:numel(subj_list)
    subj_dir = fullfile(subj_list(i).folder, subj_list(i).name);
    sessions = dir(subj_dir);
    sessions = sessions(contains({sessions.name}, 'ses'));
    for j = 1:numel(sessions)
        session_dir = fullfile(sessions(j).folder, sessions(j).name);
        anat_dir = fullfile(session_dir, 'anat');
        all_results= dir(anat_dir);
        
        % GM
        smth_str = all_results(startsWith({all_results.name}, 'smwmwp1'));
        smth = fullfile(smth_str.folder, smth_str.name);
        src = smth;
        dst_dir = fullfile(smth_dir, subj_list(i).name, sessions(j).name);
        if ~isfolder(dst_dir)
            mkdir(dst_dir)
        end
        dst = fullfile(dst_dir, smth_str.name);
        movefile(src, dst);
        
        % WM
        smth_str = all_results(startsWith({all_results.name}, 'smwmwp2'));
        smth = fullfile(smth_str.folder, smth_str.name);
        src = smth;
        dst_dir = fullfile(smth_dir, subj_list(i).name, sessions(j).name);
        if ~isfolder(dst_dir)
            mkdir(dst_dir)
        end
        dst = fullfile(dst_dir, smth_str.name);
        movefile(src, dst);
        
        % CSF
        smth_str = all_results(startsWith({all_results.name}, 'smwmwp3'));
        smth = fullfile(smth_str.folder, smth_str.name);
        src = smth;
        dst_dir = fullfile(smth_dir, subj_list(i).name, sessions(j).name);
        if ~isfolder(dst_dir)
            mkdir(dst_dir)
        end
        dst = fullfile(dst_dir, smth_str.name);
        movefile(src, dst);   
    end
end
    