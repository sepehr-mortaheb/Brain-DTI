clc
clear

%% Initialization

preproc_dir = '/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/preprocessed/ROS';
denoised_dir = '/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/denoised/ROS';

func_TR = 2;

% --- Set the Participants and Sessions Information --- 

subjectDirs = dir(fullfile(preproc_dir, 'sub-*'));
subjects = struct();
for i=1:length(subjectDirs)
    SubjName = subjectDirs(i).name;
    SubjPath = fullfile(preproc_dir, SubjName);
    
    % Get session folders
    sessionDirs = dir(fullfile(SubjPath, 'ses-*'));
    validSessions = {};

    for j=1:length(sessionDirs)
        sessName = sessionDirs(j).name;
        funcPath = fullfile(SubjPath, sessName, 'func');
        if isfolder(funcPath)
            funcContents = dir(funcPath);
            funcFiles = funcContents(~ismember({funcContents.name}, {'.', '..'}));
            if ~isempty(funcFiles)
                validSessions{end+1} = sessName;
            end
        end
    end

    subjects(i).id = SubjName;
    subjects(i).path = SubjPath;
    subjects(i).sessions = validSessions; 
end

%% Run the Denoising Loop 
for i = 1:length(subjects)
    subj = subjects(i).id;
    sessions = subjects(i).sessions;
    for j =1:length(sessions)
        ses = sessions{j};
        if isfolder(fullfile(preproc_dir, subj, ses))
            % =====================================================================
            % reading files
            
            % reading smoothed functional data in the MNI space
            fname = dir(fullfile(preproc_dir, subj, ses, ...
            'func/sra*.nii'));
            func_file = cellstr(fullfile(fname.folder, fname.name));
            % reading structural data in the MNI space
            sname = dir(fullfile(preproc_dir, subj, ses, 'anat/sub*bc.nii'));
            struct_file = cellstr(fullfile(sname.folder, sname.name));
            % reading GM mask in the MNI space
            gname = dir(fullfile(preproc_dir, subj, ses, 'anat/mri/mwp1*.nii'));
            gm_file = cellstr(fullfile(gname.folder, gname.name));
            % reading WM mask in the MNI space
            wname = dir(fullfile(preproc_dir, subj, ses, 'anat/mri/mwp2*.nii'));
            wm_file = cellstr(fullfile(wname.folder, wname.name));
            % reading CSF mask in the MNI space
            cname = dir(fullfile(preproc_dir, subj, ses, 'anat/mri/mwp3*.nii'));
            csf_file = cellstr(fullfile(cname.folder, cname.name));
            
            % reading movement regressors
            movname = dir(fullfile(preproc_dir, subj, ses, ...
                'func/rp_*.txt'));
            mov_file = cellstr(fullfile(movname.folder, movname.name));
            % reading outlier volumes
            outname = dir(fullfile(preproc_dir, subj, ses, ...
                'func/art_regression_outliers_swra*.mat'));
            out_file = cellstr(fullfile(outname.folder, outname.name));
            
            
            % =====================================================================
            % CONN batch initialization
            batch.filename = fullfile(fullfile(preproc_dir, subj, ses, 'conn_temp.mat'));
            
            % =====================================================================
            % CONN Setup
            batch.Setup.nsubjects=1;
            batch.Setup.functionals{1}{1} = func_file;
            batch.Setup.structurals{1} = struct_file;
            batch.Setup.RT = func_TR;
            batch.Setup.conditions.names = {'rest'};
            batch.Setup.conditions.onsets{1}{1}{1} = 0;
            batch.Setup.conditions.durations{1}{1}{1} = inf;
            batch.Setup.masks.Grey{1} = gm_file;
            batch.Setup.masks.White{1} = wm_file;
            batch.Setup.masks.CSF{1} = csf_file;
            batch.Setup.covariates.names = {'movement'; 'outliers'};
            batch.Setup.covariates.files = {mov_file; out_file};
            batch.Setup.analyses = 2;
            batch.Setup.analysisunits = 2;
            batch.Setup.isnew = 1;
            batch.Setup.done = 1;
            
            % =====================================================================
            % CONN Denoising
            batch.Denoising.filter=[0.008, 0.09];
            batch.Denoising.detrending = 1;
            batch.Denoising.confounds.names = {'White Matter'; 'CSF'; ...
                'movement'; 'outliers'};
            batch.Denoising.confounds.deriv = {0, 0, 1, 0};
            batch.Denoising.confounds.dimensions{1} = 5;
            batch.Denoising.confounds.dimensions{2} = 5;
            batch.Denoising.done=1;
            
            % =====================================================================
            % running batch
            conn_batch(batch);
            
            % =====================================================================
            % converting the output to the nifti format and saving it to the save
            % directory
            curr_path = pwd;
            cd(fullfile(preproc_dir, subj, ses, 'conn_temp/results/preprocessing/'));
            conn_matc2nii
            if ~isfolder(fullfile(denoised_dir, subj, ses))
                mkdir(fullfile(denoised_dir, subj, ses));
            end
            movefile('niftiDATA*.nii', fullfile(denoised_dir, subj, ses));
            delete(fullfile(preproc_dir, subj, ses, 'conn_temp.mat'));
            rmdir(fullfile(preproc_dir, subj, ses, 'conn_temp*'), 's');
            cd(curr_path)
            
            % =====================================================================
            % removing unnecessary files and directories
            
        end
    end
end