function out = func_ReadFiles(subj, ses, task_name)

ffiles = dir(fullfile(subj.dir, ses, 'func'));
ffiles = ffiles(3:end);
sfiles = dir(fullfile(subj.dir, ses, 'anat'));
sfiles = sfiles(3:end);

fdata = ffiles(contains({ffiles.name}, subj.name) & ...
               contains({ffiles.name}, ses) & ...
               contains({ffiles.name}, ['task-' task_name]) & ...
               contains({ffiles.name}, 'bold.nii'));
           
sdata = sfiles(contains({sfiles.name}, subj.name) & ...
               contains({sfiles.name}, ses) & ...
               contains({sfiles.name}, 'run-fmri') & ...
               contains({sfiles.name}, 'T1w.nii'));
           
if isempty(sdata) % if fMRI session T1 does not exist, use the DWI session T1
    sdata = sfiles(contains({sfiles.name}, subj.name) & ...
                   contains({sfiles.name}, ses) & ...
                   contains({sfiles.name}, 'run-dti') & ...
                   contains({sfiles.name}, 'T1w.nii'));
end

out = {fullfile(subj.dir, ses, 'func', fdata.name), ...
       fullfile(subj.dir, ses, 'anat', sdata.name)};