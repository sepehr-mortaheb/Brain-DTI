clc
clear 
close 

seg_dir = '/Users/sepehrmortaheb/MyDrive/LEIA/Projects/BRAIN-DTI/Astronauts-VBM/results/segmentation';
save_dir = fullfile(seg_dir, 'ready_for_python');

if ~isfolder(fullfile(seg_dir, 'ready_for_python'))
    mkdir(fullfile(seg_dir, 'ready_for_python'));
end

atlases = {'neuromorphometrics', 'Schaefer2018_100Parcels_17Networks_order'};
    
flist = dir(seg_dir);
flist = flist(contains({flist.name}, 'catROI'));

for f = 1:numel(flist)
    load(fullfile(flist(f).folder, flist(f).name));
    [~, fname, ~] = fileparts(fullfile(flist(f).folder, flist(f).name));
    for i = 1:numel(atlases)
        atlas = atlases{i};
        Rnames = S.(atlas).names;
        RVgm = S.(atlas).data.Vgm;
        RVwm = S.(atlas).data.Vwm;
        if strcmp(atlas, 'neuromorphometrics')
            RVcsf = S.(atlas).data.Vcsf; 
            save(fullfile(save_dir, [fname '_' atlas '.mat']), 'Rnames', 'RVcsf', 'RVgm', 'RVwm');
        else 
            save(fullfile(save_dir, [fname '_' atlas '.mat']), 'Rnames', 'RVgm', 'RVwm');
        end
    
    end
end





