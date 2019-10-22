% This is a matlab demo showing how to extract PET radiomic features 
% from volume of interest; this code include 57 common texture 
% features, 31 novel texture features.
% 
% Xu et al. Sub-regional radiomics analysis of PET/CT imaging with intra-tumor
% partitioning: application to prognosis for nasopharyngeal carcinoma. Molecular
% Imaging and Biology, 2019.
%
% ljlubme@gmail.com
% Southern Medical University
%
% Ackonwledgements:
% Saeed Ashrafinia: https://github.com/ashrafinia/SERA
% Martin Valli¨¨res: https://github.com/mvallieres/radiomics
%


clear all;close all;clc;
rng('default');
warning('off');

%% add the path of relative functions
addpath('...');

%% load data ( for example)
load PET_CT_data.mat

%% Parameter of Radiomics Framework Settings 
DataType  =  'PET';         % Type of the dataset. Choose from 'PET', 'CT' or 'MRscan'
pixelW  =  0.98;            % the original voxel width
sliceTh  =  3;              % the original voxel thickness

IsotVoxSize = 1;       % New isotropic voxel size for resampling in 3D. This will be the new voxel size in X, Y and Z dimension. 
VoxInterp  = 'linear';      % Image resampling interpolation type  ('nearest', 'linear', or 'cubic'). Note: 'cubic' yeilds inconsistensies with IBSI results. 
ROIInterp  = 'linear';     % ROI resampling interpolation type  ('nearest', 'linear', or 'cubic'), default: 'linear'
ROI_PV  = 0.5;          % (default 0.5) ROI partial volume threshold. Used to threshold ROI (mask) after resampling: i.e. ROI(ROI<ROI_PV) = 0, ROI(ROI>ROI_PV) = 1.

DiscType  = 'FBN';        % Discretization type: either 'FBN' (fixed bin numbers) or 'FBS' (fixed bin size or fixed bin width). 
qntzAlg  = 'Uniform';    % Discretization Type: Either 'Uniform' quantization or 'Lloyd' for Max-Lloyd quantization. (defualt: Uniform)
Nbin = 128;                   % Number of bins (for FNB) or bin size (the size of each bin for FBS). It can be an array, and the features will be calculated for each NB or BS. 

isIsot2D  = 0;            % (default 0) whether to resample image to isotropic 2D voxels (=1, i.e.keep the original slice thickness) or resample to isotropic 3D voxels (=0). (This is for 1st order features. Higher order 2D features are always calculated with original slice thickness). 
isScale  = 0;            % whether to do scaling. Has to be 1 to perform any resampling. If 0, always uses the original voxel dimension. 
isGLrounding  = 0;            % whether to round voxel intensities to the nearest integer (usually =1 for CT images, =0 for PET and SPECT)
isReSeg = 0;               % whether to perform range re-segmentation. The range is defined below in ReSegIntrvl. NOTE: Re-segmentation generally cannot be provided for arbitrary-unit modalities (MRI, SPECT)
ResegIntrval = [];   % resegmentation interval. Intensity values outside this interval would be replaced by NaN. 
isOutliers  = 0;            % whether to perform intensity outlier filtering re-segmentaion: remove outlier intensities from the intensity mask. If selected, voxels outside the range of +/- 3 standard deviation will be removed. 
isclcGlobPeak = 0;          % % whether to calculate global peak. Global peak intensity searches the whole image for peak intensity. Requires 3D resampling of the whole img and the ROI+extensive search: very time consuming. Suggest not to do for big images. If =0, then local peak = global peak. Note: in this case, global peak would be a redundant feature, should be removed from analysis.

%% initialization
RawImg = PET{1};
ROI = Mask{1};
RawImg((ROI==0)) = NaN;

%% First-order features
% Generate 3D ROI for first order and morphological features (ROI preperation)
[ImgBox, ~] = getImgBox(RawImg,ROI,isReSeg,ResegIntrval);
[SUVmax,SUVpeak,SUVmean,SUVstd,SUVvar,SUVenergy,AUC_CSH] = getSUVmetrics(ImgBox);

% ROI preperation
[ROIBox3D,levels3D,ROIonlyMorph3D,IntsROI,RawInts,RawROI,newPixW,newSliceTh]...
    = prepareVolume(RawImg,ROI,DataType,pixelW,sliceTh,1,IsotVoxSize,VoxInterp,...
    ROIInterp,ROI_PV,'XYZscale',isIsot2D,isScale,isGLrounding,DiscType,qntzAlg,...
    Nbin,isReSeg,ResegIntrval,isOutliers, isclcGlobPeak);

[texturesH] = getGlobalTextures(ROIBox3D,Nbin);
Mean_hist = texturesH.Mean;
Variance_hist = texturesH.Variance;
Skewness_hist = texturesH.Skewness;
Kurtosis_hist = texturesH.Kurtosis;
Energy_hist = texturesH.Energy;
Entropy_hist  = texturesH.Entropy;   

%% Texture Features  
[GLCM] = getGLCM(ROIBox3D,levels3D);
[texturesC] = getGLCMtextures(GLCM);

[GLRLM] = getGLRLM(ROIBox3D,levels3D);
[texturesRL] = getGLRLMtextures(GLRLM);

[GLSZM] = getGLSZM(ROIBox3D,levels3D);
[texturesSZ] = getGLSZMtextures(GLSZM);

[NGTDM,countValid] = getNGTDM(ROIBox3D,levels3D);
[texturesTD] = getNGTDMtextures(NGTDM,countValid);

[GLGLM] = getGLGLM(ROIBox3D,levels3D);
[texturesGL] = getGLGLMtextures(GLGLM);

NGLDM = getNGLDM(ROIBox3D,ROIonlyMorph3D);
texturesLD = getNGLDMtexture(NGLDM);

[TS,TS_Si] = getTS(ROIBox3D,ROIonlyMorph3D);
texturesTS = getTStexture(TS_Si);

TFC = getTFC(ROIBox3D,ROIonlyMorph3D);
texturesTFC = getTFCtexture(TFC);

[TFC,TFCCM] = getTFCCM(ROIBox3D,ROIonlyMorph3D);
texturesTFCCM = getTFCCMtexture(TFC,TFCCM);
