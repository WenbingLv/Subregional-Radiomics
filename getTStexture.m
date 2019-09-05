function textures = getTStexture(TS_Si)
%% black white symmetry
TSbws = 0;
    for idx = 1:364
        TSbws = TSbws + abs(TS_Si(idx)-TS_Si(idx+365));
    end
    TSbws = (1-TSbws/sum(TS_Si(:)))*100;    
    textures.BWS = TSbws;   

%% Max spectrum
    Maxspe = max(TS_Si(:))/sum(TS_Si(:));
    textures.MaxSpe = Maxspe;
