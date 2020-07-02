function CombiningEstimates(ImputationsDirectory)
% Script for combining estimates for multiple imputations. The output is a
% folder called "MeanAnalysis" in which results can be interrogated using
% SPM, and contrasts specified. "MultipleImputations_SPMest" needs to be run first, then second level
% analysis run on the imputed data sets.  Inputs to this script is the
% directory in which these second level analyses are held.  Folders for
% each analysis should begin with the prefix "Imp"
%
% Original by Elizabeth Liddle & Marjie Jansen; modified by Paul Briley on 22/3/2019

mkdir(ImputationsDirectory,'MeanAnalysis');
DestinationDirectory = [ImputationsDirectory filesep 'MeanAnalysis'];

ImputedFolders = dir([ImputationsDirectory filesep 'Imp*']);
Nimputations = size(ImputedFolders,1);

MaskDir = dir([ImputationsDirectory filesep ImputedFolders(1).name filesep 'mask*']);
[~,FILENAME,EXT]=fileparts(MaskDir(1).name);
if ~strcmp(EXT,'.nii'); EXT='.img'; end

V = spm_vol([ImputationsDirectory filesep ImputedFolders(1).name filesep FILENAME EXT]);
Y = spm_read_vols(V);
V.fname = [DestinationDirectory filesep FILENAME EXT];
V.private.dat.fname = V.fname;
spm_write_vol(V,Y);

RPVDir = dir([ImputationsDirectory filesep ImputedFolders(1).name filesep 'RPV*']);
[~,FILENAME,EXT]=fileparts(RPVDir(1).name);
if ~strcmp(EXT,'.nii'); EXT='.img'; end

V = spm_vol([ImputationsDirectory filesep ImputedFolders(1).name filesep FILENAME EXT]);
Y = spm_read_vols(V);
V.fname = [DestinationDirectory filesep FILENAME EXT];
V.private.dat.fname = V.fname;
spm_write_vol(V,Y);

% Assembles info about the imputed beta images.

BetaDir = dir([ImputationsDirectory filesep ImputedFolders(1).name filesep 'beta*' EXT]);
Nbetas = size(BetaDir,1);

CombinedImageVar = zeros(size(Y));

% Computes and writes mean beta images for each regressor, and computes the
% between imputation variance.  Also computes mean resels and
% FWHM parameters and writes them into the new SPM.mat file
for iBeta = 1:Nbetas
    Images = zeros(V.dim(1),V.dim(2),V.dim(3),Nimputations);
    for iImputation = 1:Nimputations
        [~,BetaName,EXT] = fileparts(BetaDir(iBeta).name);
        if ~strcmp(EXT,'.nii'); EXT = '.img'; end
        
        ImpBetaPath = [ImputationsDirectory filesep ImputedFolders(iImputation).name filesep BetaName EXT];        
        V = spm_vol(ImpBetaPath);
        Y = spm_read_vols(V);
        Images(:,:,:, iImputation) = Y;
    end
    CombinedImageMean = mean(Images,4);
    
    RepmatMeans = zeros(size(Images));    
    for iImputations = 1:Nimputations   
        RepmatMeans(:,:,:,iImputations) = CombinedImageMean;
    end
    CombinedImageVar = sum((Images - RepmatMeans).^2,4)./(Nimputations-1);
    
    Vbeta = V;
    Vbeta.fname = [DestinationDirectory filesep BetaName EXT];
    Vbeta.private.dat.fname = Vbeta.fname;
    
    spm_write_vol(Vbeta,CombinedImageMean)
end

% Computes the within imputation variance and writes a new ResMS image
% representing combined error variance
ResMSs = zeros(V.dim(1),V.dim(2),V.dim(3),Nimputations);
for iImputation = 1:Nimputations    
    ResMSImagePath = [ImputationsDirectory filesep ImputedFolders(iImputation).name filesep 'ResMS' EXT];
    
    VresMS = spm_vol(ResMSImagePath);
    YresMS = spm_read_vols(VresMS);
    ResMSs(:,:,:,iImputation) = YresMS;
    
    load([ImputationsDirectory filesep ImputedFolders(iImputation).name filesep 'SPM.mat'],'SPM');
    R(iImputation,:) = SPM.xVol.R;
    FWHM(iImputation,:) = SPM.xVol.FWHM;
    
    clear SPM;
end

MeanR = mean(R);
MeanFWHM = mean(FWHM);

OldSPMPath = [ImputationsDirectory filesep ImputedFolders(1).name filesep 'SPM.mat'];
NewSPMPath = [DestinationDirectory filesep 'SPM.mat'];

copyfile(OldSPMPath,NewSPMPath);

load(NewSPMPath,'SPM');

SPM.xVol.R = MeanR;
SPM.xVol.FWHM = MeanFWHM;

nSubj = length(SPM.xY.P);
for s = 1:nSubj
    for iImputation = 1:Nimputations
        [PATHSTR,FILENAME,EXT] = fileparts(SPM.xY.P{s}(1:end-2));
        % FILENAME in format: ImpNo1_con_XXXX
        if ~strcmp(FILENAME(1:6),'ImpNo1'); error('filename not in expected format'); end
        FILENAME(6) = num2str(iImputation);
        
        S = spm_vol([PATHSTR filesep FILENAME EXT]);
        T = spm_read_vols(S);
        if iImputation==1
            myMeanImage = T;
        else
            myMeanImage = myMeanImage + T;
        end
    end
    myMeanImage = myMeanImage./Nimputations;
    
    pos = nan;
    for z = 1:length(FILENAME)
        if isnan(pos) && strcmp(FILENAME(z),'_')
            pos = z;
            S.fname = [PATHSTR filesep 'ImpMean' FILENAME(pos:end) EXT];
            S.private.dat.fname = S.fname;
        end
    end
    if isnan(pos); error('filename not as expected'); end
    
    spm_write_vol(S,myMeanImage);
    
    SPM.xY.VY(s).fname = S.fname;
    SPM.xY.VY(s).private.dat.fname = S.fname;
    SPM.xY.P{s} = [S.fname ',1'];
end    

save(NewSPMPath,'SPM');

MeanResMS = mean(ResMSs,4);

ResDir = dir([ImputationsDirectory filesep ImputedFolders(1).name filesep 'ResMS*']);
[~,FILENAME,EXT] = fileparts(ResDir(1).name);
TotalVariance = (1+1./Nimputations) .* CombinedImageVar + MeanResMS;
NewResMSPath = [DestinationDirectory filesep FILENAME EXT];

VnewResMS = VresMS;
VnewResMS.fname = NewResMSPath;
VnewResMS.private.dat.fname = NewResMSPath;
spm_write_vol(VnewResMS,TotalVariance)

% Computes the contrasts for the combined analysis
temp = SPM.xCon;
consize = length(temp);
SPM.swd = DestinationDirectory;
cd(DestinationDirectory)
save SPM SPM

for iCon = 1:consize
    spmtemp = spm_FcUtil('Set', temp(iCon).name, temp(iCon).STAT, 'c', temp(iCon).c, SPM.xX.xKXs);
    SPM.xCon(iCon) = spmtemp;
end
if consize > 0
    spm_contrasts(SPM);
end

disp('Your new analysis is now ready in the folder "MeanAnalysis"');
