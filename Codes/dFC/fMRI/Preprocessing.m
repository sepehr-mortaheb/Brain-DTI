clc 
clear 

%% Initialization 

% --- Set the following directories --- 

% Directory of the BIDS formated data:
bids_dir = '/Users/sepehrmortaheb/Desktop/pipeline_test/data';
% Save directory of the fMRI processing:
save_dir = '/Users/sepehrmortaheb/Desktop/pipeline_test/preprocessed';

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

% Subjects list [Ex: {'sub-XXX'; 'sub-YYY'}]
subj_list = {'sub-control01'};

% Sessions list [Ex: {'ses-ZZZ'; 'ses-TTT'}]
ses_list = {'ses-01'};

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

% Subject Information Struct
Subjects(length(subj_list)) = struct();
for i=1:length(subj_list)
    Subjects(i).name = subj_list{i};
    Subjects(i).dir = fullfile(bids_dir, subj_list{i});
    Subjects(i).sessions = ses_list; 
end

% Adding required paths 
addpath(art_dir);
addpath(spm_dir);
addpath(fullfile(spm_dir, 'src'));
addpath('./functions');

%% Functional Pipeline 

for subj_num = 1:numel(subj_list)
    subj = subj_list{subj_num};
    func_PipelineSS(Dirs, Subjects(subj_num), AcqParams);
end