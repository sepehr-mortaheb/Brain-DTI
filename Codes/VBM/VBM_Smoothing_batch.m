function matlabbatch = VBM_Smoothing_batch(data)

matlabbatch = cell(1,1);
matlabbatch{1}.spm.spatial.smooth.data = data;
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';