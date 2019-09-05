function textures = getTFCCMtexture(TFC,TFCCM)
 %%   
vec = TFC(:);
vec = vec(find(vec>0));

mean1 = mean(vec);
sd1 = std(vec);

tempvar = 0;
    for idx1 = 1:max(TFC(:))
        for idx2 = 1:max(TFC(:))
            if TFCCM(idx1,idx2)>0
                tempvar = tempvar  - TFCCM(idx1,idx2)*log(TFCCM(idx1,idx2))  ;
            end
        end
    end
    textures.CodeEntropy = tempvar;

    %%
    tempvar = 0;
    for idx1 = 1:max(TFC(:))
        %for idx2 = 1:max(FC3D_global(:))
            if TFCCM(idx1,idx1)>0
                tempvar = tempvar + TFCCM(idx1,idx1)^2  ;
            end
        %end
    end
    textures.CodeSimilarity = tempvar;

%%
    if length(find(isnan(TFCCM))) > 1
        feature_output = 0;
    else
        feature_output = 0;
        for idx1 = 1:size(TFCCM,1)
            for idx2 = 1:size(TFCCM,2)
                feature_output = feature_output + TFCCM(idx1, idx2)*(idx1-idx2)^2;
            end
        end
    end
textures.Contrast = feature_output;

%%
textures.SAM = sum(sum(TFCCM.^2)); % The Second angular moment is a summation of squared co-occurrence matrix

%%
    tempvar = graycoprops1(round(TFCCM),'Correlation');
    textures.Correlation = tempvar.Correlation;
 %%
     feature_output = 0;
    for idx1 = 1:size(TFCCM,1)
        for idx2 = 1:size(TFCCM,2)
            feature_output = feature_output + TFCCM(idx1, idx2)*(1/(1+(idx1-idx2)^2));
        end
    end
    textures.IDM = feature_output;
 
 %%
     feature_output = 0;
    if length(find(isnan(TFCCM))) < 1     
        for idx1 = 1:size(TFCCM,1)
            for idx2 = 1:size(TFCCM,2)
                feature_output = feature_output + TFCCM(idx1, idx2)*(1/(1+abs(idx1-idx2)));
            end
        end
    end
    textures.Homogeneity = feature_output;

    %%
    feature_output = 0;
    for idx1 = 1:size(TFCCM,1)
        for idx2 = 1:size(TFCCM,2)
            feature_output = feature_output + idx1*idx2*TFCCM(idx1, idx2);
        end
    end
    textures.Intensity = feature_output;
    
    %%
    textures.Entropy = -sum(sum(TFCCM.*log(TFCCM + realmin))); 

    