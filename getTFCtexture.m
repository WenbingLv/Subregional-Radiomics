function textures = getTFCtexture(TFC)
        
    textures.Coarseness = length(find(TFC==20))/ length(find(TFC>0)) ;
    
    textures.Homogeneity = length(find(TFC==1))/ length(find(TFC>0)) ;
%%    
vec = TFC(:);
vec = vec(find(vec>0));

mean1 = mean(vec);
sd1 = std(vec);

    tempvar = 0;
    for idx1 = 1:max(TFC(:))
        tempvar = tempvar + abs(idx1*length(find(TFC == idx1))/numel(vec)-mean1);
    end
    textures.MeanCovergence = tempvar/sd1 ;
    
  %%     
    tempvar = 0;
    for idx1 = 1:max(TFC(:))
        tempvar = tempvar + (idx1-mean1)^2*length(find(TFC == idx1))/numel(vec);
    end
    textures.Variance = tempvar;
    
    %%


    
