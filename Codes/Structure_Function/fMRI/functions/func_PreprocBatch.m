function matlabbatch = func_PreprocBatch(inpfiles, AcqParams, Dirs)

% Setting some parameters 
stc_num = AcqParams.nslc;
stc_ord = AcqParams.ordslc;
stc_ref = AcqParams.refslc;
tr = AcqParams.tr;

spm_dir = Dirs.spm;
fdata = inpfiles{1};
sdata = inpfiles{2};

% Defining the slice order variable 
switch stc_ord
    case 1
        slice_order = [1:1:stc_num];
    case 2
        slice_order = [stc_num:-1:1];
    case 3
        for k=1:stc_num
            slice_order = round((stc_num-k)/2 + (rem((stc_num-k), 2) * (stc_num-1)/2)) + 1;
        end
    case 4
        slice_order = [1:2:stc_num 2:2:stc_num];
    case 5
        slice_order = [stc_num:-2:1 stc_num-1:-2:1];
    case 6
        fname = fdata;
        fname(end-2:end)=[];
        fname = [fname 'json'];
        fid = fopen(fname); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        slice_order = val.SliceTiming; 
end

%% Reading the structural and functional data 
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'struct';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {cellstr(sdata)}';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'func';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {cellstr(fdata)};


%% Slice Timing Correction
matlabbatch{3}.spm.temporal.st.scans{1}(1) = cfg_dep('Named File Selector: func(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{3}.spm.temporal.st.nslices = stc_num;
matlabbatch{3}.spm.temporal.st.tr = tr;
matlabbatch{3}.spm.temporal.st.ta = tr - (tr/stc_num);
matlabbatch{3}.spm.temporal.st.so = slice_order;
matlabbatch{3}.spm.temporal.st.refslice = stc_ref;
matlabbatch{3}.spm.temporal.st.prefix = 'a';

%% Realignment
matlabbatch{4}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%% Coregistration of Functional Data to the T1 Space 
matlabbatch{5}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Named File Selector: struct(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{5}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{5}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

%% Normalization of Functional Data to the MNI Space 
% matlabbatch{7}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
% matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
% matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-Inf -Inf -Inf
%                                                           Inf Inf Inf];
% matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
% matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;
% matlabbatch{7}.spm.spatial.normalise.write.woptions.prefix = 'w';

%% Smoothing
matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{6}.spm.spatial.smooth.dtype = 0;
matlabbatch{6}.spm.spatial.smooth.im = 0;
matlabbatch{6}.spm.spatial.smooth.prefix = 's';
