function func_PipelineSS(Dirs, subj, AcqParams)

out_dir = Dirs.out;
task_name = AcqParams.name;
ses_list = subj.sessions;

%% Initialization
if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% Preprocessing loop over the sessions 
for i = 1:numel(ses_list)
    ses = ses_list{i};
    
    % If the data exist in this session run the pipeline
    if length(dir(fullfile(subj.dir, ses, 'anat'))) > 2
        % reading the data
        inpfiles = func_ReadFiles(subj, ses, task_name);
    
        % perform preprocessing
        func_Preproc(inpfiles, Dirs, subj, ses, AcqParams)
        close all
    end
end
