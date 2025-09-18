clc 
clear 

mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI';
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results/SC';

atlas = 'SCH100';
R = 100;

evoi = [1, 2, 3, R-2, R-1, R];

mni_node = readtable(fullfile(mni_dir, [atlas '.csv']));
load(fullfile(res_dir, ['SC_eig_' atlas '.mat']));

mat = zeros(R, 5);

for r=evoi
    mat(:, 1) = mni_node.R;
    mat(:, 2) = mni_node.A;
    mat(:, 3) = mni_node.S;
    mat(:, 4) = U(:, r);
    mat(:, 5) = ones(R, 1);
    writematrix(mat, fullfile(res_dir, ['SC_eig_' num2str(r) '_' atlas '.txt']),...
        'Delimiter', 'tab');
end
    


