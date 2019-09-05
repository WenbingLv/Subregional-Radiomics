function textures = getNGLDMtexture(NGLDM)

%% Entropy
NGLDM_entropy = 0;
    for idx1 = 1:size(NGLDM,1)
        for idx2 = 1:size(NGLDM,2)
            if NGLDM(idx1,idx2)>0
                NGLDM_entropy = NGLDM_entropy + NGLDM(idx1,idx2)*log(NGLDM(idx1,idx2));
            end
        end
    end
    textures.Entropy = - NGLDM_entropy / sum(NGLDM(:));

%% Second moment(Energy)
    NGLDM_Energy = sum(sum(NGLDM.^2));    
    textures.Energy = NGLDM_Energy / sum(NGLDM(:));

%% Small number emphasis
NGLDM_SNE = 0;
    for idx1 = 1:size(NGLDM,1)
        for idx2 = 1:size(NGLDM,2)
            NGLDM_SNE = NGLDM_SNE + NGLDM(idx1, idx2)/idx2^2;
        end
    end
    textures.SNE = NGLDM_SNE / sum(NGLDM(:));
    
%% Large number emphasis
NGLDM_LNE = 0;
    for idx1 = 1:size(NGLDM,1)
        for idx2 = 1:size(NGLDM,2)
            NGLDM_LNE = NGLDM_LNE + NGLDM(idx1, idx2) * idx2^2;
        end
    end
    textures.LNE = NGLDM_LNE / sum(NGLDM(:));
    
%% number nonuniformity
NGLDM_NNU = 0;
    for idx2 = 1:size(NGLDM,2)        
        NGLDM_NNU = NGLDM_NNU + sum(NGLDM(:, idx2))^2;
    end
    textures.NNU = NGLDM_NNU / sum(NGLDM(:));
