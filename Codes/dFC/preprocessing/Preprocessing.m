clc 
clear 

%% Initialization 

% --- Set the following directories --- 

% Directory of the BIDS formated data:
bids_dir = '/Volumes/LACIE_SHARE/BRAIN-DTI_dFC/data/bids_ros';
% Save directory of the fMRI processing:
save_dir = '/Volumes/LACIE_SHARE/BRAIN-DTI_dFC/data/preprocessed';

%##########################################################################
% --- Set the Acquisition Parameters --- 

% ------ General Parameters ------
% The name of the functional task
task_name = 'rest';
% Repetition Time (RT) of the functional acquisition (seconds)
func_TR = 2; 

% ------ STC-Specific Parameters ------ 
% Number of Slices
stc_num = 42;
% Slice Order (1=ascending, 2=descending, 3=interleaved(middle-top),
% 4=interleaved(buttom-up), 5=interleaved(top-down), 6:slice timings
% available in the JSON file)
stc_ord = 6;
% Reference Slice/Time
stc_ref = 0;

%##########################################################################
% --- Set the Participants Information --- 

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
        funcPath = fullfile(SubjPath, sessName, 'func');
        if isfolder(funcPath)
            funcContents = dir(funcPath);
            funcFiles = funcContents(~ismember({funcContents.name}, {'.', '..'}));
            if ~isempty(funcFiles)
                validSessions{end+1} = sessName;
            end
        end
    end

    Subjects(i).name = SubjName;
    Subjects(i).dir = SubjPath;
    Subjects(i).sessions = validSessions; 
end

%##########################################################################
% --- Creating Handy Variables and AddPath Required Directories ---

% Directories Struct
art_dir = which('art');
art_dir(end-4:end) = []; 
spm_dir = which('spm');
spm_dir(end-4:end) = [];
Dirs = struct();
Dirs.bids = bids_dir; 
Dirs.out = save_dir;
Dirs.spm = spm_dir;
Dirs.art = art_dir;

% Acquisition Parameters Struct
AcqParams = struct();
AcqParams.name = task_name;
AcqParams.tr = func_TR; 
AcqParams.nslc = stc_num;
AcqParams.ordslc = stc_ord;
AcqParams.refslc = stc_ref;

% Adding required paths 
addpath(fullfile(spm_dir, 'src'));
addpath('./functions');

%% Functional Pipeline 

for subj_num = 1:numel(Subjects)
    func_PipelineSS(Dirs, Subjects(subj_num), AcqParams);
end