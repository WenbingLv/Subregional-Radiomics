% this is an example of how to split whole tumor into several sub-regions;
% this code mainly contained three steps: local entropy maps generation, 
% individual-level clustering, population-level clustering.
%
% ljlubme@gmail.com
% Southern Medical University
%
clear all;close all;clc;
rng('default');
warning('off')
%% load data
load PET_CT_data.mat
%% add path
addpath('...')
%% parameter setting
K_cluster_1 = 40;
patch = ones(9,9,9);
%% local entropy maps generation
for i=1:length(PET)
    % generate PET/CT local entropy maps
    PET_box = getSUVbox(PET{i},Mask{i});
    CT_box = getSUVbox(CT{i},Mask{i});
    Mask_box = getSUVbox(Mask{i},Mask{i});
    PET_box = normalization(PET_box);
    CT_box = normalization(CT_box);
    PET_entropy_box = entropyfilt(PET_box,patch);
    CT_entropy_box = entropyfilt(CT_box,patch);
    % set out-of-tumor pixels are NaNs
    CT_box(~Mask_box)=NaN;
    PET_box(~Mask_box)=NaN;
    PET_entropy_box(~Mask_box) = NaN;
    CT_entropy_box(~Mask_box) = NaN;
    % normalization to keep four paramaters with balanced contibution 
    PET_entropy_box = normalization(PET_entropy_box);
    CT_entropy_box = normalization(CT_entropy_box);
    % incorprated into 4D vector
    data_all{i} = [PET_box(:),CT_box(:),PET_entropy_box(:),CT_entropy_box(:)]; 
    % save all maps
    PET_box_all{i} = PET_box;
    CT_box_all{i} = CT_box;
    PET_entropy_all{i} = PET_entropy_box;
    CT_entropy_all{i} = CT_entropy_box;
    Mask_box_all{i} = Mask_box;
end
%% individual-level cluster
superpixel_all = [];
for cases = 1:length(data_all)
    data = data_all{cases};%for patient i
    [IDX,C] = kmeans(data,K_cluster_1,'Distance','sqEuclidean');
    IDX_all_1{cases} = IDX; % cluster label
    for k = 1:K_cluster_1
        idi = find(IDX==k);
        sp = data(idi,:);
        mean_sp = nanmean(sp,1);
        superpixel(k,:) = mean_sp;
    end
    superpixel_all = [superpixel_all;superpixel];
end
%% population-level cluster
superpixel_all = zscore(superpixel_all);
eva = evalclusters(superpixel_all,'Kmeans','CalinskiHarabasz','KList',[2:10]);
figure,plot(eva);
IDX_all_2 = mat2cell(eva.OptimalY,ones(length(PET),1).*40,1);
for cases = 1:length(PET)
    superpixels_label = IDX_all_2{cases};% superpixels label
    pixels_label = IDX_all_1{cases};% pixels label
    IDX_temp = pixels_label;
    for i=1:K_cluster_1
        IDX_temp(pixels_label==i) = superpixels_label(i);%Consistently labeled supervoxels within each tumor were merged as a sub-region  
    end
    subregion{cases} = reshape(IDX_temp,size(PET_box_all{cases}));
end
%% save results
for i=1:size(subregion{1},3)
    figure,imshow(subregion{1}(:,:,i),[]);colormap jet;caxis([0 eva.OptimalK]);
end
save('subregion.mat','subregion');
save('Maps_box_all.mat','PET_box_all','CT_box_all','PET_entropy_all','CT_entropy_all','Mask_box_all');
