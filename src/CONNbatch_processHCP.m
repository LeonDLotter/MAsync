% Script to process HCP data for INS MA ROI-to-ROI RSFC analyses

%% SETTINGS: 
PROJECTname='conn_HCP.mat';
TARGETpath='/media/leon/data_m2/MAsync/RSFC/conn'; % target folder for conn project 
CONNECTOMEpath='/media/leon/Data_4TB/RSFC/HCP_data_resampled/%s'; % HCP data
MASKpath='/media/leon/data_m2/MAsync/RSFC/vols/brainmask_mni152_3mm.nii'; % analysis space
ROIpath='/media/leon/data_m2/MAsync/RSFC/vols/macm_rTPJ_idx_forCONN.nii'; % macm ROIs
ROIname='macm';
RUNPARALLEL=true; % run parallel
NSUBJECTS=[]; % number of subjects to include, empty for all
NJOBS=4; % number of parallel jobs (empty for one job per subject)
COPYFILES=true;  % copy files to TARGETpath?
OVERWRITE=false;  % overwrites files if they exist in target folder 
overwr=0; % overwrite files in each processing steps


%% FINDS STRUCTURAL/FUNCTIONAL FILES
clear FUNCTIONAL_FILE* STRUCTURAL_FILE;
subs=dir(regexprep(CONNECTOMEpath,'%s.*$','*')); 
subs=subs([subs.isdir]>0);
subs={subs.name};
subs=subs(cellfun(@(s)all(s>='0'&s<='9'),subs));
if isempty(NSUBJECTS), NSUBJECTS=numel(subs); 
else subs=subs(1:NSUBJECTS);
end
if isempty(NJOBS), NJOBS=NSUBJECTS; end
NJOBS=min(NSUBJECTS,NJOBS);

for n=1:numel(subs)
    fprintf('Locating subject %s files\n',subs{n});
    
    % anatomical
    t1=fullfile(sprintf(CONNECTOMEpath,subs{n}),'T1','T1w_restore_brain.nii'); % T1  
    % resting-state, FIX-denoised
    f1=fullfile(sprintf(CONNECTOMEpath,subs{n}),'rest','rfMRI_REST1_LR_hp2000_clean_3mm.nii'); 
    f2=fullfile(sprintf(CONNECTOMEpath,subs{n}),'rest','rfMRI_REST1_RL_hp2000_clean_3mm.nii'); 
    if isempty(dir(t1)), error('file %s not found',t1); end
    if isempty(dir(f1)), error('file %s not found',f1); end
    if isempty(dir(f2)), error('file %s not found',f2); end
    
    % copy files if "true"
    if COPYFILES
        fprintf('Copying files to local folder\n');
        [ok,nill]=mkdir(TARGETpath,'LocalCopyDataFiles');
        [ok,nill]=mkdir(fullfile(TARGETpath,'LocalCopyDataFiles'),subs{n});
        t1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'structural.nii');  if OVERWRITE||isempty(dir(t1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',t1,t1b)); end; t1=t1b;
        f1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional1.nii'); if OVERWRITE||isempty(dir(f1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f1,f1b)); end; f1=f1b;
        f2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional2.nii'); if OVERWRITE||isempty(dir(f2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f2,f2b)); end; f2=f2b;
    end
    
    STRUCTURAL_FILE{n,1}=t1;
    FUNCTIONAL_FILE(n,1:2)={f1,f2};
end
nsessions=2;
fprintf('%d subjects, %d sessions\n',NSUBJECTS,nsessions);

% Project mat file path
projectMat = fullfile(TARGETpath,PROJECTname);


%% CREATES CONN BATCH STRUCTURE
clear batch;
batch.filename=projectMat;
if RUNPARALLEL
    batch.parallel.N=NJOBS; % number of parallel processing batch jobs, default profile
end                         

%% CONN Setup                                           
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS; % nr of subj
batch.Setup.RT=0.72; % TR

% conditions
batch.Setup.conditions.names={'rest'}; % single condition (aggregate across all sessions)
for ncond=1
    for nsub=1:NSUBJECTS
        for nses=1:nsessions   
            batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; 
            batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;
        end
    end
end     % rest condition (all sessions)

% set functional files
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]); 
for nsub=1:NSUBJECTS
    for nses=1:nsessions
        batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE(nsub,nses);
    end
end

% set structural
batch.Setup.structurals=STRUCTURAL_FILE;

% voxel size
batch.Setup.voxelmask=1; % only for voxel-wise analyses, explicit mask
batch.Setup.voxelmaskfile=MASKpath; % mni152 3mm brainmask
batch.Setup.voxelresolution=1; % same as mask & functionals (3mm isotropic)

% ROIs
batch.Setup.rois.names={ROIname};
batch.Setup.rois.files={ROIpath};
batch.Setup.rois.dimensions={1}; % average time series
batch.Setup.rois.weighted=0; % not weight by mask value
batch.Setup.rois.multiplelabels=1; % associated text file
batch.Setup.rois.dataset=0; % compute from "primary" dataset
batch.Setup.rois.add=0; % add this to the existing gm/csf/wm rois

% final
batch.Setup.analyses=[1,2]; % ROI-to-ROI & Seed-to-voxel
batch.Setup.overwrite=overwr;                            
batch.Setup.done=1;

% run batch
conn_batch(batch);

%% CONN Denoising     
clear batch;
batch.Denoising.filter=[0.01, 0.08]; % band-pass filter
batch.Denoising.detrending=1; % linear detrending
batch.Denoising.despiking=0; % no despiking
batch.Denoising.confounds.names={''}; % no additional confound regression
batch.Denoising.overwrite=overwr;
batch.Denoising.done=1; 

% run batch
conn_batch(batch);

%% CONN Analysis   
clear batch;
batch.Analysis.name={'RtR'};
batch.Analysis.sources={[ROIname '.']}; % all macm rois
batch.Analysis.measure=2; % semipartial correlations
batch.Analysis.weight=1;
batch.Analysis.modulation=0;
batch.Analysis.type=1;
batch.Analysis.overwrite=overwr;
batch.Analysis.done=1;

% run batch
conn_batch(batch);
