function [] = SourceSpaceFCanalysisBang_EloretaOnly_ICAEdited(cfg, participantnumber);
% This program was copoed from the SourceSpaceFCanalysisBang_EloretaOnly.m and works only on the ICA Edited Data;
%% Define parameters and methods;                   
if ~isfield(cfg, 'foi'),              cfg.foi              = 'alpha';                 end
if ~isfield(cfg, 'fcmethod'),         cfg.fcmethod         = 'wpli';                  end
if ~isfield(cfg, 'experimentnumber'), cfg.experimentnumber = 2;                       end
if ~isfield(cfg, 'atlastype'),        cfg.atlastype        = 'LPBA';                  end 
if ~isfield(cfg, 'replace'),          cfg.replace          = 0;                       end
if ~isfield(cfg, 'parcmethod'),       cfg.parcmethod       = 'average';               end
if ~isfield(cfg, 'methodtype'),       cfg.methodtype       = 'eloreta';               end
if ~isfield(cfg, 'mmtype'),           cfg.mmtype           = '5mm';                   end
if ~isfield(cfg, 'hztype'),           cfg.hztype           = 1;                       end
if ~isfield(cfg, 'permutation')       cfg.permutation      = 1;                       end
if ~isfield(cfg, 'skullcond'),        cfg.skullcond        = 0.0132;                  end
foi = cfg.foi; fcmethod = cfg.fcmethod; experimentnumber = cfg.experimentnumber; atlastype=cfg.atlastype; replace = cfg.replace; 
parcmethod = cfg.parcmethod;  methodtype = cfg.methodtype; mmtype = cfg.mmtype; hztype=cfg.hztype;permutation = cfg.permutation;skullcond = cfg.skullcond;

%% Global variables;
global EEG_FT EEG_freq source sourcedata moments m roitrialdata roitrialdata1 source_conn 
global elec filter atlas roivol
ft_defaults;
%if ~exist('CURRENTERP','var');eeglab;close;end;

%% Load the time series data and define MRI number;
%% Load the subjectsmatrix
if ismac
    load('/Users/wanzexie/Documents/MATLAB/WX/MATLAB/data/LCN/BangladeshProject/NSR/subjectsmatrix.mat','subjectsmatrix');
else
    load('C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\NSR\subjectsmatrix.mat','subjectsmatrix');
end
participantnumbers = subjectsmatrix.participantnumbers;
age                = subjectsmatrix.ages(participantnumbers == participantnumber);
ntrial             = subjectsmatrix.goodntrials(participantnumbers == participantnumber);
%% Check the num of trials and the existence of the output file;
if ntrial < 60 
    disp(['Participant ' num2str(participantnumber) ' has less than 60 trials']);
    return
end
if age > 12 
    segmentedtype='segmented'; %'segmented refers to the FEM model;
else
    segmentedtype='segmented_nma';%'segmentednma refers to the FEM model but with non-myelinated axons;
end
switch age   
    case 6
        mrinumber = 302; 
    case 36
        mrinumber = 328;
end
if mrinumber>=1000
    mristring=num2str(mrinumber);
elseif mrinumber>100
    mristring = ['0' num2str(mrinumber)];
else 
    mristring =  ['00' num2str(mrinumber)];
end
%Frequency of interest;
bpfreq = foi;

% check if the output already exists;
if ismac
    filepath = '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/data/LCN/BangladeshProject/Fieldtrip/fcanalysis_source/';
    sourcepath  = '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/data/LCN/BangladeshProject/Fieldtrip/sourcedata/'
    modelspath  = '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/models/';
    atlasespath = '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/atlases/';
    elecfolder = '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/Electrodes/';
    roidatafilepath ='/Users/wanzexie/Documents/MATLAB/WX/MATLAB/data/LCN/BangladeshProject/Fieldtrip/roidata/';
    eegfilepath =  '/Users/wanzexie/Documents/MATLAB/WX/MATLAB/data/LCN/BangladeshProject/Fieldtrip/rawdata/';
else
    filepath = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\Fieldtrip\fcanalysis_source\';
    sourcepath  = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\Fieldtrip\sourcedata\';
    modelspath  = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\models\';
    atlasespath = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\atlases\';
    elecfolder  = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\Electrodes\';
    roidatafilepath ='C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\Fieldtrip\roidata\';
    eegfilepath = 'C:\Users\ch197255\Documents\MATLAB\WX\MATLAB\data\LCN\BangladeshProject\Fieldtrip\rawdata\';
end
if skullcond == 0.0132
    outputname = [filepath 'Participant ' num2str(participantnumber) ' fcanalysis_' methodtype '_' num2str(foi(1)) 'to' num2str(foi(2)) 'Hz_' parcmethod '_' mmtype '_' fcmethod '_' atlastype '_' num2str(hztype) 'Hz_ica1.mat'];
else
    outputname = [filepath 'Participant ' num2str(participantnumber) ' fcanalysis_' methodtype '_' num2str(foi(1)) 'to' num2str(foi(2)) 'Hz_' parcmethod '_' mmtype '_' fcmethod '_' atlastype '_' num2str(skullcond*1000) 'cond_' num2str(hztype) 'Hz_ica1.mat'];
end
% check line95 to modify the name of the outputname;
if exist(outputname) & replace == 0;
    disp(['The following file already exists:' outputname])
    return
end
%check if the ROI data exists already; It means nROI x time series (e.g.,56 x 500 (amples) for the LPBA atlas);
if skullcond == 0.0132
    roidataname = [roidatafilepath 'Participant ' num2str(participantnumber) ' fcanalysis_' methodtype '_' parcmethod '_' mmtype '_' atlastype '_' num2str(hztype) 'Hz_ica1' '_ROIs.mat'];
else
    roidataname = [roidatafilepath 'Participant ' num2str(participantnumber) ' fcanalysis_' methodtype '_' parcmethod '_' mmtype '_' atlastype '_' num2str(skullcond*1000) 'cond_' num2str(hztype) 'Hz_ica1' '_ROIs.mat'];
end
if exist(roidataname);
    load(roidataname,'roitrialdata');
    switch atlastype
        case 'LPBA'
            nroi = 56;
        case 'IXI'
            nroi = 83;
    end
else %if no ROI based data, the do the FC analysis from the raw datasets, which means source analysis needs to be done;
%load ICA Edited EEG data; 
%For source analysis, use the non-overlapping data and maybe create overlapped trials later;
eegfilename = ['Experiment ' num2str(experimentnumber) ' Subject ' num2str(participantnumber) ' Fieldtrip ' num2str(hztype) 'Hz Highpass_ICA_Edited.mat']; 
load ([eegfilepath eegfilename],'EEG_FT'); 
%change the FT format from preprocessing to timelock format, which makes some of the following analyses easier;
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
EEG_FT = ft_timelockanalysis(cfg,EEG_FT);

%% Cortical source reconstruction;
mrifieldtripfoldername =  modelspath;
% If use eLORETA then there is no need to do source analysis from the scratch. Instead, use the filter created for each average MRI;
% However, if use beamforming then the source analysis needs to be conducted for each dataset because
% beamforming source analysis calculates the covariance matrix for each dataset.;
% The name of the electrodetype needs to be defined regardless of the methodtype for source analysis;
electrodetype='hydrocelgsn124';
%load the grid; %The grid will be used later for atlas stuff regardless of source analysis methods;
filename=['S' mristring '_sourcemodel_' mmtype '_ft_grid.mat'];
if ~exist([mrifieldtripfoldername filename]);
 disp(['does not exist quick return ' mrifieldtripfoldername filename]);
 return
end
 disp(['load ' mrifieldtripfoldername filename]);
 load([mrifieldtripfoldername filename],'grid');
 
filedirectory_lf = mrifieldtripfoldername;


if skullcond == 0.0132;
    filename_filter=['S' mristring '_' electrodetype '_' segmentedtype '_' mmtype '_ft_' methodtype '_filter.mat'];
else
    filename_filter=['S' mristring '_' electrodetype '_' segmentedtype '_' mmtype '_ft_' num2str(skullcond*1000) 'cond_' methodtype '_filter.mat'];
end
%S0302_hydrocelgsn124_segmented_nma_5mm_ft_66cond_eloreta_filter

if exist([filedirectory_lf filename_filter]);
    disp(['load ' filedirectory_lf filename_filter]);
    load([filedirectory_lf filename_filter],'filter');
    efilter = filter;
else
    disp(['no saved filter quick return ' filedirectory_lf filename_filter]);
    return
end

%Then do source reconstruction for each single trial using the filter created in the previous source analysis;
sourcedata = [];
ntrial  = size(EEG_FT.trial,1);
ndipole = find(~cellfun(@isempty,efilter)); %This array equals find(grid.inside);
%ndipole = find(grid.inside);
ntime   = size(EEG_FT.time,2);
sourcedata.trial = NaN(ntrial,length(ndipole),ntime);
tic;
for i=1:ntrial
    fprintf('%.0f out of %.0f \n',i,ntrial); %Inform the process (n out of the total trial number);
    for j = 1:length(ndipole)
        moments = efilter{ndipole(j)} *squeeze(EEG_FT.trial(i,:,:));        
        [u, s, v] = svd(moments, 'econ');  %Find the dimension explaining/representing the most variance; 
        m = u(:,1)'*moments; %this m equals to or has a linear relationship with v(:,1)', e.g., all the elements are 3 times bigger than the elements in v(:,1)';
        %Or simply do:  %m = v(:,1)'; 
        sourcedata.trial(i,j,:) = m;    
    end
end
analysistime = toc;
disp(['The Source analysis took ' num2str(analysistime) 's']);
%save([sourcepath sourcefilename],'sourcedata');

%% Do PCA or calculate the average CDR for each brain region;
tic;
%load the atlas;
switch atlastype
    case 'LPBA'
        atlasname = ['S' mristring '_LPBA40atlas_ft.mat'];
        disp(['load ' atlasespath atlasname]);
        load([atlasespath atlasname],'atlas');
        nroi = 56;
    case 'IXI'
        disp(['load ' atlasespath ixiANTSatlas_ft]);
        load([atlasespath ixiANTSatlas_ft],'atlas');
        nroi = 83;
end

%load the atlas/brain ac info
brain_ac = ['S' mristring '_brain_atlas_ac.mat'];
load([modelspath brain_ac],'brain_ac');
%Find the dipoles with activity;
positioninside=grid.pos(find(grid.inside),:);
positioninside=positioninside+double(int64(repmat(brain_ac,size(positioninside,1),1)));
indicesinside=sub2ind(atlas.dim,positioninside(:,1),positioninside(:,2),positioninside(:,3));

%Do the PCA for each brain ROI;
for j=1:length(atlas.labelnumbers)
 searchloc=find(atlas.inside == atlas.labelnumbers(j));    %indices for this label
 searchpos=atlas.pos(searchloc,:);                 %locations for this label
 searchpos=double(int64(searchpos))+repmat(brain_ac,size(searchpos,1),1);    
 searchindices=sub2ind(atlas.dim,searchpos(:,1),searchpos(:,2),searchpos(:,3));
 if strcmp(parcmethod,'centroid')
     % find out the centroid position;
     centroidpos = round(mean(searchpos));
     % create a centroid cube;
     X = [centroidpos(1)-10:centroidpos(1)+10];
     Y = [centroidpos(2)-10:centroidpos(2)+10];
     Z = [centroidpos(3)-10:centroidpos(3)+10];
     cubepos = combvec(X,Y,Z)';
     cubeindices=sub2ind(atlas.dim,cubepos(:,1),cubepos(:,2),cubepos(:,3));
     %get an array that marks every point in sourcedata with atlas (cubes)
     [jj kk]=ismember(indicesinside,cubeindices);
     [junk1 junk2]=ismember(indicesinside,searchindices);
     junk3 = sum([jj junk1],2);
     junk3 == 2;
     jj(junk3~=2) = 0;
 else
      %get an array that marks every point in sourcedata with atlas (ROIs)
     [jj kk]=ismember(indicesinside,searchindices);
 end
 atlasmatrix(j,:)=jj';  
end
roivol=atlasmatrix*repmat(1,[size(atlasmatrix,2) 1]);    %sum over columns
atlasmatrixdiv=atlasmatrix./repmat(roivol,[1 size(atlasmatrix,2)]); % This is for the average method;

roitrialdata=[];
roitrialdata.trial = NaN(ntrial,length(atlas.labelnumbers),ntime);
for i = 1:ntrial;
    switch parcmethod
        case 'PCA'
            for j = 1:length(atlas.labelnumbers)
                roidata = [];
                roidata = sourcedata.trial(i,find(atlasmatrix(j,:)),:);
                roidata = squeeze(roidata)';
                [coeff score latent useless explained] = pca(roidata);
                roidata_maxpca =  score(:,1)'; % pick up the pca component explaining the most of the variance;
                %{ 
                roidata_pcas = roidata*coeff;  %original component, and the scores are these components minus the means;
                roidata_maxpca2 = roidata_pcas(:,1)'; 
                %} 
                roitrialdata.trial(i,j,:)=roidata_maxpca;
            end
        case {'average','centroid'}
            roitrialdata.trial(i,:,:) = atlasmatrixdiv*squeeze(sourcedata.trial(i,:,:));
    end
end

roitrialdata.time = EEG_FT.time;
roitrialdata.dimord = 'rpt_chan_time';
for k = 1:length(ndipole);
    roitrialdata.label = atlas.labels;
end;
if exist('sourcedata','var'); clear sourcedata;end
analysistime = toc;
disp(['The averaging/PCA/centroid process took ' num2str(analysistime) 's']);

%Save the roitrialdata;
disp(['Save: ' roidataname]);
save(roidataname,'roitrialdata');

end %if exist(roidatname);

%% Do freqanalysis and FC analysis
EEG_freq = []; source_conn = []; EEG_FT = roitrialdata;
%do permutation to control of observation (trial num) bias for wpli;
if permutation & strcmp(fcmethod,'wpli')
lowernumber   = 80;
biggernumber  = ntrial;
a = zeros(50,lowernumber);
for i = 1:50
    a(i,:) = sort(randperm(biggernumber,lowernumber));
end
fcmatrix  = zeros(size(a,1),nroi,nroi); 
fcmatrixz = zeros(size(a,1),nroi,nroi); 
for i = 1:size(a,1);
EEG_FT_pm = EEG_FT;
fprintf('%.0f out of %.0f \n',i,size(a,1)); %inform process
    EEG_FT_pm.trial = EEG_FT_pm.trial(a(i,:),:,:); % this EEG_FT format is preprocessing not nonprocessing fieldtrip format;
    EEG_FT_pm.time  = EEG_FT_pm.time;
    % freq analysis
    cfg           = [];
    cfg.method    = 'mtmfft';
    if min(bpfreq)<6
    cfg.pad       = 1.024; %http://www.bitweenie.com/listings/fft-zero-padding/;
    end
    cfg.taper     = 'hanning'; %cfg.taper = 'dpss'; 'dpss' simply means apply multiple tapers;
    cfg.output    = 'fourier';
    cfg.foilim    = bpfreq;
    EEG_freq = ft_freqanalysis(cfg, EEG_FT_pm);
    %take the average of the fourier coefficients over the frequencies for an foi (e.g., 6, 7, 8, 9hz for alpha; Miskovic & Keil, 2015);           
    %EEG_freq.fourierspctrm = mean(EEG_freq.fourierspctrm,3); %calculate the mean fourier coefficients;
    %EEG_freq.freq = mean(cfg.foilim);
    %EEG_freq.dimord = 'rpttap_chan';        

    % functional connectivity analysis;
    cfg           = [];
    cfg.method    = fcmethod;
    source_conn = ft_connectivityanalysis(cfg, EEG_freq);
    %If the diagonal values are 1 then change them to 0. This needs to be checked for different methods; 
    %  if strcmp(fcmethod,'coh');
    %      source_conn.cohspctrm = source_conn.cohspctrm - eye(length(source_conn.cohspctrm));
    %  end
    %DO fisher z transform; This step is highly recommended when compare the coherence/correlation values statistally.
    switch fcmethod
        case 'wpli'
            source_conn.wplispctrm = abs(source_conn.wplispctrm); %use the abs value for wpli;
            
            source_conn.wplispctrm = mean(source_conn.wplispctrm,3);
            
            source_conn.zwplispctrm = 0.5*log((1+source_conn.wplispctrm)./(1-source_conn.wplispctrm));
            fcmatrix(i,:,:) = source_conn.wplispctrm;
            fcmatrixz(i,:,:)= source_conn.zwplispctrm;
            %figure;imagesc(source_conn.wplispctrm); 
        case 'wpli_debiased'
            source_conn.zwpli_debiasedspctrm = 0.5*log((1+source_conn.wpli_debiasedspctrm)./(1-source_conn.wpli_debiasedspctrm));
            fcmatrix(i,:,:) = source_conn.wpli_debiasedspctrm;
            fcmatrixz(i,:,:)= source_conn.zwpli_debiasedspctrm;
            %figure;imagesc(source_conn.wpli_debiasedspctrm); 
    end
end 
fcmatrix  = squeeze(mean(fcmatrix,1));
fcmatrixz = squeeze(mean(fcmatrixz,1));

else % no permutation for the other methods;
%freq analysis
cfg            = [];
cfg.output     = 'fourier';
%cfg.method = {'coh', 'csd', 'wpli', 'wpli_debiased', 'plv','imag'}
cfg.method     = 'mtmfft';
if min(bpfreq)<6
cfg.pad        = 1.024;
end
cfg.foilim     = bpfreq;
cfg.keeptrials = 'yes';
cfg.taper      = 'hanning';
EEG_freq       = ft_freqanalysis(cfg, EEG_FT);
%take the average of the fourier coefficients over the frequencies for an foi (e.g., 6, 7, 8, 9hz for alpha; Miskovic & Keil, 2015);           
%EEG_freq.fourierspctrm = mean(EEG_freq.fourierspctrm,3); %calculate the mean fourier coefficients;
%EEG_freq.freq = mean(cfg.foilim);
%EEG_freq.dimord = 'rpttap_chan';      

%fc analysis
cfg = [];
cfg.method = fcmethod;
if strcmp(fcmethod,'imag');
    cfg.method = 'coh';
    cfg.complex = 'absimag';
end
source_conn = ft_connectivityanalysis(cfg, EEG_freq);
switch fcmethod
    case {'coh','imag'}
        source_conn.cohspctrm = mean(source_conn.cohspctrm,3);
        source_conn.zcohspctrm = 0.5*log((1+source_conn.cohspctrm)./(1-source_conn.cohspctrm));
        fcmatrix(:,:) = source_conn.cohspctrm;%calculate the mean fc for different frequency bins;
        fcmatrixz(:,:)= source_conn.zcohspctrm;
end

end % if permutation & strcmp(fcmethod,'wpli')
% put the fcmatrix and fcmatrix into the source_conn structure;
fcmatrix(logical(eye(size(fcmatrix))))   = 0;
fcmatrixz(logical(eye(size(fcmatrixz)))) = 0;
%figure;imagesc(fcmatrix); imagesc(fcmatrixz); 
source_conn = [];
source_conn.fcmatrix  = fcmatrix;
source_conn.fcmatrixz = fcmatrixz;

%% Save the source_conn matrix;
disp(['Save: ' outputname]);
save(outputname,'source_conn');
