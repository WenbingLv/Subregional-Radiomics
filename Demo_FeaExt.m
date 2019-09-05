%This is a matlab demo showing how to extract radiomic features 
% from volume of interest; this code include 57 common texture 
% features, 31 novel texture features.
%
% ljlubme@gmail.com
% Southern Medical University
%
clear all;close all;clc;
rng('default');
warning('off');

%% add path
addpath('...\Texture-Function');

%% load data
load Maps_box_all.mat

%% parameter setting
pixelsize = 0.98; 
sliceS = 3;
gray_level = 128;

ROIonlyH = PET_box_all{1};
ROIonlyS = Mask_box_all{1};
              
%% First-order features
[SUVmax,SUVpeak,SUVmean,SUVstd,SUVvar,SUVenergy,AUC_CSH] = getSUVmetrics(ROIonlyH);

[texturesH] = getGlobalTextures(ROIonlyH,128);
Mean_hist = texturesH.Mean;
Variance_hist = texturesH.Variance;
Skewness_hist = texturesH.Skewness;
Kurtosis_hist = texturesH.Kurtosis;
Energy_hist = texturesH.Energy;
Entropy_hist  = texturesH.Entropy;   

%% Texture Features
[ROIonlyM,levelsM] = uniformQuantization(ROIonlyH,gray_level);
   
[GLCM] = getGLCM(ROIonlyM,levelsM);
[texturesC] = getGLCMtextures(GLCM);

[GLRLM] = getGLRLM(ROIonlyM,levelsM);
[texturesRL] = getGLRLMtextures(GLRLM);

[GLSZM] = getGLSZM(ROIonlyM,levelsM);
[texturesSZ] = getGLSZMtextures(GLSZM);

[NGTDM,countValid] = getNGTDM(ROIonlyM,levelsM);
[texturesTD] = getNGTDMtextures(NGTDM,countValid);

[GLGLM] = getGLGLM(ROIonlyM,levelsM);
[texturesGL] = getGLGLMtextures(GLGLM);

NGLDM = getNGLDM(ROIonlyM,ROIonlyS);
texturesLD = getNGLDMtexture(NGLDM);

[TS,TS_Si] = getTS(ROIonlyM,ROIonlyS);
texturesTS = getTStexture(TS_Si);

TFC = getTFC(ROIonlyM,ROIonlyS);
texturesTFC = getTFCtexture(TFC);

[TFC,TFCCM] = getTFCCM(ROIonlyM,ROIonlyS);
texturesTFCCM = getTFCCMtexture(TFC,TFCCM);
