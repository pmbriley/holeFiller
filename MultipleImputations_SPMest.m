function MultipleImputations_SPMest(OriginalDirectory,Threshold,NumberImputations)
% Run an initial 2nd level analysis on your subjects' con images. Then
% enter the directory of this analysis
%
% USAGE:
% MultipleImputations_SPMest(OriginalDirectory,Threshold,NumberImputations)
%
% OriginalDirectory = directory containing the original second level analysis
% Threshold = proportion of NaNs allowed per voxel (i.e. proportion of missing subjects allowed), typically use 0.15
% NumberImputations = number of imputations (e.g. 5)
%
% Original by Elizabeth Liddle & Marjie Jansen; modified by Paul Briley on 22/3/2019

spm('defaults','fmri');

ImputationsDirectory = [OriginalDirectory '_imputed'];
mkdir(ImputationsDirectory);

% load the SPM.mat file from the original 2nd level analysis
load(fullfile(OriginalDirectory,'SPM.mat'),'SPM');
OriginalSPM = SPM; clear SPM;

% X is the design matrix from the original analysis
X = OriginalSPM.xX.X;
FileCells = OriginalSPM.xY.P;
FileNames = char(FileCells);
if strcmp(FileNames(1,end-1:end),',1')
    FileNames = FileNames(:,1:end-2);
end

pKX = OriginalSPM.xX.pKX;
xKXs = OriginalSPM.xX.xKXs;
Nbetas = size(OriginalSPM.xX.Bcov,2);
Vc = zeros(Nbetas,1);

Bmat = zeros(Nbetas,1);
for iBeta= 1:Nbetas
    Bmat(iBeta,1) = 1;
    Vc(iBeta,1) = Bmat'*OriginalSPM.xX.Bcov*Bmat;
end

% For a between-subjects analysis, this will be the number of subjects
Nimages = size(FileNames,1);

if Threshold<1
    Threshold=floor(Nimages.*Threshold);
end

W = speye(Nimages,Nimages);

nFilledHoles = zeros(1,NumberImputations);
imputedVoxels = cell(1,NumberImputations);
for iImputation = 1:NumberImputations
    disp('busy with imputation number ');
    disp(iImputation);
    
    ImpDirectory = [ImputationsDirectory filesep 'ImpNo' num2str(iImputation)];
    mkdir(ImpDirectory);
    
    for iImage = 1:Nimages % subjects' con images        
        [ImagePathString,ImageName,EXT]=fileparts(FileNames(iImage,:));
        ImagePath=[ImagePathString filesep ImageName EXT];        
        if iImage==1           
            V = spm_vol(ImagePath);            
            Y=spm_read_vols(V);
            Dim1 = V.dim(1);
            Dim2 = V.dim(2);
            Dim3 = V.dim(3);
            Images = zeros(Nimages, Dim1, Dim2, Dim3);
            Vnew = V;           
        else            
            V = spm_vol(ImagePath);
            Y = spm_read_vols(V);
        end        
        Images(iImage,:,:,:) = Y;
    end    

    NewImages = Images; % Images is a 4-D matrix of 3D vols with subjects the 1st Dim
    
    for iDim1 = 1:Dim1
        for iDim2 = 1:Dim2
            for iDim3 = 1:Dim3                
                y = Images(:,iDim1,iDim2,iDim3); % find values in same voxel across subject
                Holes = find(isnan(y)); % find where the holes are
                NHoles = numel(Holes); % number of holes
                
                if NHoles<=Threshold && NHoles>0 % if voxel is to be filled
                    nFilledHoles(iImputation) = nFilledHoles(iImputation) + 1;
                    fillVec = [iDim1; iDim2; iDim3];
                    imputedVoxels{iImputation} = [imputedVoxels{iImputation} fillVec];                    
                    
                    GoodIndex = isnan(y)==false;
                    Good_y = y(GoodIndex);
                    NGoodValues = numel(Good_y);
                    
                    Wvoxel      = W(GoodIndex,GoodIndex);
                    
                    Good_pKX    = pKX(:,GoodIndex);
                    KWY         = spm_filter(OriginalSPM.xX.K,Wvoxel*Good_y);
                    xX.xKXs     = OriginalSPM.xX.xKXs;
                    xX.xKXs.X   = xX.xKXs.X(GoodIndex,:);
                    xX.xKXs.u   = xX.xKXs.u(GoodIndex,:);
                    
                    b           = Good_pKX*KWY; % Parameter estimates
                    res         = spm_sp('res',xX.xKXs)*KWY; % Residuals           
                    
                    ResSS       = sum(res.^2);
                    SE          = sqrt(ResSS*Vc);
                    MSE         = ResSS/(NGoodValues-length(b));
                    ErrorVariance=sqrt(MSE);
                    
                    % Residual SSQ
                    clear KWY
                    
                    RandBetas = nan(1,length(b));
                    for m = 1:length(b)
                        RandBetas(1,m) = normrnd(b(m),SE(m));
                    end
                    
                    % matrix with Beta coefficients for columns and number of holes for rows
                    RandBetasMat   = repmat(RandBetas,NHoles,1);
                    HolePredictors = X(isnan(y),:);
                    HolePredicteds = sum(HolePredictors.*RandBetasMat,2);
                    
                    Imputations = nan(1,NHoles);
                    for n = 1:NHoles
                        Imputations(1,n) = HolePredicteds(n)+normrnd(0,ErrorVariance);
                    end
                    
                    NewImages(isnan(y),iDim1,iDim2,iDim3) = Imputations;
                    clear Imputations;                   
                end
            end
        end
    end
    
    SPM = OriginalSPM;    
    
    P = cell(Nimages,1);
    for n = 1:Nimages       
        [PATHSTR,NAME,EXT]=fileparts(FileNames(n,:));           
            
        if iImputation==1
            mkdir([PATHSTR '_imputed']);
        end
        
        NewImagePathName = [PATHSTR '_imputed' filesep 'ImpNo' num2str(iImputation) '_' NAME EXT];        
        
        Vnew.fname = NewImagePathName;
        Vnew.private.dat.fname = NewImagePathName;
        I = squeeze(NewImages(n,:,:,:));
       
        spm_write_vol(Vnew,I);
        
        P{n,1} = [NewImagePathName ',1'];
        
        SPM.xY.VY(n).fname = Vnew.fname;
        SPM.xY.VY(n).private.dat.fname = Vnew.fname;        
    end
    
    MeanImage = squeeze(mean(NewImages,1));
    MaskImage = NaN(size(MeanImage));
    MaskImage(~isnan(MeanImage)) = 1;
    
    V.fname = [ImpDirectory filesep 'MeanImage.nii'];
    V.private.dat.fname = V.fname;
    spm_write_vol(V,MeanImage);
    
    V.fname = [ImpDirectory filesep 'xmask.nii'];
    V.private.dat.fname = V.fname;
    spm_write_vol(V,MaskImage);    
    
    cd(ImpDirectory);
 
    SPM.xY.P = P;
    cd(ImpDirectory);
    SPM.swd = ImpDirectory;
    save SPM SPM
    temp = SPM.xCon;
    consize = length(SPM.xCon);
   
    SPM = rmfield(SPM,{'Vbeta','VResMS','VM','xCon','SPMid','xVol'});
    SPM.xX = rmfield(SPM.xX, {'K','W', 'xKXs','pKX', 'V', 'trRV', 'trRVRV', 'erdf', 'Bcov','nKX'});
    
    toRemove = {'UFp','h','Cy'};
    for rem = 1:length(toRemove)
        try SPM.xVi = rmfield(SPM.xVi,toRemove{rem});
        catch; end
    end
    
    SPM = spm_spm(SPM);
    SPM.xCon = temp;
    
    for j = 1:consize
        spmtemp = spm_FcUtil('Set', temp(j).name, temp(j).STAT, 'c', temp(j).c, SPM.xX.xKXs);
        SPM.xCon(j) = spmtemp;
    end
    
    if consize > 0
        spm_contrasts(SPM);
    end
end % Imputation loop

save(fullfile(ImputationsDirectory,'FilledHoles.mat'),'nFilledHoles','imputedVoxels');
CombiningEstimates(ImputationsDirectory);
