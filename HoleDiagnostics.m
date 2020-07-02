function SmallMasks = HoleDiagnostics
% This function analyses the con images you entered into your second level
% analysis for holes.  It draws a series of plots that will help you decide
% on how to approach the problem of missing data, and creates and saves a 
% series of diagnostic images that you can inspect later. It will also, if you
% want, output the filenames of any subjects with concerningly small masks.
% If you want to output these as a variable, use the syntax
% VarName=HoleDiagnostics; either way, it will display them in the command
% window

% The first graphic to be output is an SPM display of a brain image
% weighted voxel wise by the  proportion of subjects for whom data is
% availabe for that voxel.  This  image can be surfed to explore the
% brain regions where most good data is available

% The next is a plot with three graphs

% The top is a bar chart, with one bar for each of your subjects, plotted
% in the order in which they were entered into your analysis.  The vertical
% axis represents the number of voxels in each subject's con image as a
% proportion of the number of voxels final collective mask (the mask called mask.img in your
% second level analysis, which represents the voxels in which all subjects
% have good data).  No subject can therefore have less than 100% of this
% number; however many subjects may have many more.  If this plot reveals
% that one or a few subjects's values are close to one, but the rest are very much larger,
% this suggests that it may be worth looking at these subjects' data, to
% ensure there isn't a fixable reason (e.g. bad realignment) that is
% making their image so small.  However if most subject's values are
% substantially more than one, it suggests that the missing values are well
% distributed throughout the subjects, and that repairing the analysis for
% these voxels may be worthwhile

% The middle plot is a raster plot in which the voxels are ranked along the
% X axis according to the proportion of subjects for which there are data,
% while on the Y axis, subjects are ranked by the number of voxels for
% which there are data.  From this plot, it should be possible to evaluate
% how evenly distributed the missing data is throughout the subjects

% The bottom plot is a sensitivity analysis, showing the potential increase
% in mask size (increase in number of voxels used in the analysis) will
% result from an increasing number of allowable NaNs for each voxel (i.e an
% increasing number of missing data points for each voxel)

% The function then offers you the opportunity to see the effect of
% allowing a given number of missing datapoints on the number of voxels
% available for analysis.  The resulting images are shown in an SPM checkreg 
% display in which the top left image shows the proportion data for each voxel; 
% the top right shows the extent of the mask after thresholding; the bottom
% left shows the original mask volume, and the bottom right shows the mean
% con image - the mean of all available con image data for each voxel.
% These images are in register and can therefore be surfed to check that
% additional brain areas that will be included are likely to be within the
% brain, and have sensible data

% Finally, you will be asked whether you want to see the filenames for any
% subjects with concerningly small masks, as indicated on the bar chart.
% If so, enter a vector of numbers from the bar chart

% When you run the function you will be prompted for the path to your
% original second level analysis

% Original by Elizabeth Liddle and Marjie Jansen; modified by Paul Briley
% on 22/3/2019

spm('defaults','fmri');
scrsz = get(0,'ScreenSize');

% This is the directory containing the original second level analysis
OriginalDirectory = uigetdir('*.*', 'Select the directory of the analysis you want to diagnose');

%This loads the SPM.mat file from the original 2nd level analysis
load([OriginalDirectory filesep 'SPM.mat'],'SPM');
X=SPM.xX.X; %X is the design matrix from the original analysis

FileNames=char(SPM.xY.P);
FileNames=FileNames(:,1:end-2);
NmaskVoxels=size(SPM.xVol.XYZ,2);
clear SPM

%For a between-subjects analysis, this will be the number of subjects
Nimages=size(FileNames,1);

for iImage=1:Nimages % subjects' con images    
    [ImagePathString,ImageName,EXT]=fileparts(FileNames(iImage,:));
    
    if strcmp(EXT, '.nii')==true % added to allow for both .nii file and Analyze files
        HeaderSuffix ='.nii';
        ImageSuffix ='.nii';
    else
        HeaderSuffix ='.hdr';
        ImageSuffix ='.img';
    end
    
    ImagePath=[ImagePathString filesep ImageName HeaderSuffix];
    
    if iImage==1
        V=spm_vol(ImagePath);
        Y=spm_read_vols(V);

        Images=zeros(size(Y,1), size(Y,2), size(Y,3),Nimages);
        Nvoxels=numel(Y);
    end
    
    V=spm_vol(ImagePath);
    Y=spm_read_vols(V);

    Images(:,:,:,iImage)=Y;
end
clear Filenames;
cd(OriginalDirectory);

% create a new image that is calculated by taking the mean of the con
% images
MeanConFilename=['NanMeanCon' ImageSuffix];
MeanConImage=nanmean(Images,4);
V=spm_vol(['ResMS' HeaderSuffix]);

Y=spm_read_vols(V);
V.fname=MeanConFilename;
V.descrip='Mean of first level con images';
spm_write_vol(V,MeanConImage);
clear MeanConImage

% create a new image that consists of proportion of available subjects per voxel
% (i.e. an image with Nans for nans and ones for non-nans)
DiagnosticFilename  = ['DiagnosticImage' ImageSuffix];
V                   = spm_vol(['ResMS' HeaderSuffix]);
V.fname             = DiagnosticFilename;
V.descrip           = 'Diagnostics: proportion of voxels present';
VoxelsPresent       = zeros(size(Images));
VoxelsPresent(isfinite(Images))=1;

nVoxelsPresent      = nansum(VoxelsPresent,4); % number of non nans per voxel
MaxPresentVoxels=max(max(max(nVoxelsPresent)));
Y   = nVoxelsPresent/MaxPresentVoxels; % as prop of total number of subjects
PropVoxelsPresent=reshape(Y,Nvoxels,1);

spm_write_vol(V,Y);
clear Images nVoxelsPresent

% Displays brain image weighted by proportion of subjects in whom each
% voxel is present
matlabbatch{1}.spm.util.disp.data = ...
    {[OriginalDirectory filesep DiagnosticFilename ',1']};

spm_jobman('run', matlabbatch)

h = gcf;
movegui(h,'northeast')

clear matlabbatch

% Draws bar chart by subject number (L-R) of voxels for each subject
figure
h=gcf;
set(h,'OuterPosition',[1 1 scrsz(3)/2 scrsz(4)])

VoxelsPerSubject=NaN(Nimages,1);
for iImage=1:Nimages
    PresentVox=find(VoxelsPresent(:,:,:,iImage)==1);   
    VoxelsPerSubject(iImage,1)=length(PresentVox);
end

PropVoxelsPerSubject=VoxelsPerSubject/MaxPresentVoxels;

subplot(3,1,1)
bar(PropVoxelsPerSubject)
axis([0 Nimages 0 max(PropVoxelsPerSubject)])
set(gca,'XTick', 1:Nimages,'XTickLabel', 1:Nimages, 'fontsize',6, 'XTickLabelRotation',90);

xlabel('Subjects');
ylabel('Number of present voxels/current mask size');

% Draws map of voxels present, ranked vertically by number of voxels in
% each subject (smallest to largest images) and horizontally by subjects
% per voxel (least to most).

VoxelsPresentReshape=reshape(VoxelsPresent,Nvoxels,Nimages);
AtLeastOneIndex=PropVoxelsPresent>0; % vector
AtLeastOneImages=VoxelsPresentReshape(AtLeastOneIndex,:); % take 1/NaN info from those with at least one voxel present, matrix
RankedVoxels=PropVoxelsPresent(AtLeastOneIndex); % take propVoxelsPresent for only >0 values, vector
AtLeastOneImages=[RankedVoxels,AtLeastOneImages]; % make a new matrix, first column propVoxelsPresent
AtLeastOneImages=sortrows(AtLeastOneImages,1); % sort the matrix by propVoxelsPresent
AtLeastOneImages=AtLeastOneImages(:,2:end)'; % take original matrix back
AtLeastOneImages=[PropVoxelsPerSubject,AtLeastOneImages]; % adds column of prop of Voxels per subject
AtLeastOneImages=sortrows(AtLeastOneImages,1); % sorts by prop voxels per subject
AtLeastOneImages=AtLeastOneImages(:,2:end)'; % take orig matrix back

subplot(3,1,2)
colormap('gray')
image(AtLeastOneImages', 'CDataMapping', 'scaled')
caxis([0, 1])
xlabel('Voxels ranked by proportion of data present');
ylabel('Subjects ranked by total number of present voxels');

% Draws another picture to show how many voxels would be included at
% certain thresholds of number of nans allowed
RemainingVoxels=zeros(Nimages,2);
for i=1:Nimages
    Threshold=i/Nimages;
    ThresholdedVoxelsPresentIndex=PropVoxelsPresent>=Threshold;
    RemainingVoxels(i,1)=Nimages-i; % threshold
    RemainingVoxels(i,2)=sum(ThresholdedVoxelsPresentIndex)./NmaskVoxels; % number of voxels included compared to mask
end

subplot(3,1,3)
scatter(RemainingVoxels(:,1),RemainingVoxels(:,2))
axis([0 max(RemainingVoxels(:,1)) 1 max(RemainingVoxels(:,2))]);
xlabel('Number of NaNs allowed');
ylabel('increase in mask size');

% Allows experimenter to see the effects on mask size of different
% thresholds for allowable missing data (.05 - allows voxels with up to 5%
% missing data, etc)
TryThresholds='y';
Counter=0;
while TryThresholds== 'y'
    Counter=Counter+1;
    
    if Counter<2
        Button=questdlg('Do you want to see the effect of a threshold value for included voxels?');
    else
        Button=questdlg('Do you want to see the effect of another threshold value for included voxels?');
    end
    
    switch Button
        case 'Yes'
            TryThresholds='y';
        case 'No'
            TryThresholds='n';
        case 'Cancel'
            TryThresholds='y';
    end
    
    if  TryThresholds=='y'
        NmissingAllowed=inputdlg('Enter the number of allowable missing subjects:');
        Threshold=str2num(char(NmissingAllowed))/Nimages;
        
        ThresholdedVoxelsPresent=PropVoxelsPresent;
        ThresholdedVoxelsPresentIndex=PropVoxelsPresent<=1-Threshold;
        ThresholdedVoxelsPresent(ThresholdedVoxelsPresentIndex)=0;
        TestDiagnosticFilename=['DiagnosticImage_' char(NmissingAllowed) '_filled.img'];
        
        V=spm_vol(['ResMS' HeaderSuffix]);
        V.fname=TestDiagnosticFilename;
        V.descrip='Diagnostics: proportion of voxels present over threshold';
        Y=reshape(ThresholdedVoxelsPresent,V.dim);
        spm_write_vol(V,Y);
        
        clear ThresholdedVoxelsPresent
        
        Cell={[OriginalDirectory filesep DiagnosticFilename ',1'];...
            [OriginalDirectory filesep TestDiagnosticFilename ',1'];...
            [OriginalDirectory filesep 'mask' ImageSuffix ',1'];...
            [OriginalDirectory filesep MeanConFilename]};
        
        matlabbatch{1}.spm.util.checkreg.data = Cell;
        spm_jobman('run', matlabbatch)
        text(-450,800,'Proportion data for each voxel')
        text(-50,800,['Thresholded at ' num2str(round(Threshold*100)) '%'])
        text(-400,100,'Original mask volume')
        text(-70,100,'mean con image from all subjects')
        
        shg
        clear matlabbatch
    end
end

Button=questdlg('Do you want to see the filenames of subjects with small masks?');
switch Button
    case 'Yes'
        SmallMaskIndex=str2num(char(inputdlg('Enter a vector of numbers from the bar chart')));
        if nargout>0
            SmallMasks=FileNames(SmallMaskIndex,:);
        else
            disp(FileNames(SmallMaskIndex,:));
        end
    case 'No'
        if nargout>0
            SmallMasks='';
        end
    case 'Cancel'
        if nargout>0
            SmallMasks='';
        end        
end
