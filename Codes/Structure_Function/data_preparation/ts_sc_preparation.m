clc
clear

%% Initialization

preproc_dir = '/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/preprocessed/ROS';
denoised_dir = '/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/denoised/ROS';
out_dir = '/home/smortaheb/Projects/Brain-DTI/Structure-Function/Data/TS_SC/ROS';

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
    subjects(i).sessions = validSessions; 
end

%% Reading TS and SC data and move them to the output directory

% Read and copy SC data 

for i = 1:length(subjects)
    subj = subjects(i).id;
    sessions = subjects(i).sessions;
    for j =1:length(sessions)
        ses = sessions{j};
        if ~isfolder(fullfile(out_dir, subj, ses))
            mkdir(fullfile(out_dir, subj, ses));
        end
        % AAL
        sc_aal_src = fullfile(preproc_dir, subj, ses, 'dwi/SC_AAL.csv');
        sc_aal_dst = fullfile(out_dir, subj, ses, 'SC_AAL.csv');
        copyfile(sc_aal_src, sc_aal_dst);

        % SCH100
        sc_aal_src = fullfile(preproc_dir, subj, ses, 'dwi/SC_SCH100.csv');
        sc_aal_dst = fullfile(out_dir, subj, ses, 'SC_SCH100.csv');
        copyfile(sc_aal_src, sc_aal_dst);

        % SCH400
        sc_aal_src = fullfile(preproc_dir, subj, ses, 'dwi/SC_SCH400.csv');
        sc_aal_dst = fullfile(out_dir, subj, ses, 'SC_SCH400.csv');
        copyfile(sc_aal_src, sc_aal_dst);
    end
end

% Read and Copy TS data 
for i = 1:length(subjects)
    subj = subjects(i).id;
    sessions = subjects(i).sessions;
    for j =1:length(sessions)
        ses = sessions{j};
        if ~isfolder(fullfile(out_dir, subj, ses))
            mkdir(fullfile(out_dir, subj, ses));
        end
        ts_mat = fullfile(denoised_dir, subj, ses, 'ROI_Subject001_Condition000.mat');
        load(ts_mat);
        data = cell2mat(data(4:end));
        ts_AAL = data(:, 501:end);
        ts_SCH100 = data(:, 1:100);
        ts_SCH400 = data(:, 101:500);

        save(fullfile(out_dir, subj, ses, 'ts_AAL.mat'), 'ts_AAL');
        save(fullfile(out_dir, subj, ses, 'ts_SCH100.mat'), 'ts_SCH100');
        save(fullfile(out_dir, subj, ses, 'ts_SCH400.mat'), 'ts_SCH400');
    end
end