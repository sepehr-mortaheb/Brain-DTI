function data = VBM_read_data(mode, data_dir, seg_dir)

if mode == 1
    copyfile(data_dir, seg_dir);
    subj_list = dir(seg_dir);
    subj_list = subj_list(contains({subj_list.name}, 'sub'));
    data = cell(numel(subj_list), 1);
    for i = 1:numel(subj_list)
        subj_dir = fullfile(subj_list(i).folder, subj_list(i).name);
        sessions = dir(subj_dir);
        sessions = sessions(contains({sessions.name}, 'ses'));
        images = cell(numel(sessions), 1);
        for j = 1:numel(sessions)
            session_dir = fullfile(sessions(j).folder, sessions(j).name);
            anat_dir = fullfile(session_dir, 'anat');
            T1s = dir(anat_dir);
            T1s = T1s(startsWith({T1s.name}, 'sub') & contains({T1s.name}, 'T1w.nii') & contains({T1s.name}, 'acq-dti'));
            if numel(T1s) > 1
                T1s = T1s(1);
            end
            images{j} = fullfile(T1s.folder, [T1s.name ',1']);
        end
        data{i} = images;
    end
    
elseif mode ==3
    subj = input('participnat name = ', 's');
    subj_dir = fullfile(data_dir, subj);
    dest_dir = fullfile(seg_dir, subj);
    if ~isfolder(dest_dir)
        mkdir(dest_dir)
    end
    copyfile(subj_dir, dest_dir);
    sessions = dir(dest_dir);
    sessions = sessions(contains({sessions.name}, 'ses'));
    images = cell(numel(sessions), 1);
    data = cell(1, 1);
    for j = 1:numel(sessions)
        session_dir = fullfile(sessions(j).folder, sessions(j).name);
        anat_dir = fullfile(session_dir, 'anat');
        T1s = dir(anat_dir);
        T1s = T1s(startsWith({T1s.name}, 'sub') & contains({T1s.name}, 'T1w.nii') & contains({T1s.name}, 'acq-dti'));
        if numel(T1s) > 1
            T1s = T1s(1);
        end
        images{j} = fullfile(T1s.folder, [T1s.name ',1']);
    end
    data{1} = images;
end