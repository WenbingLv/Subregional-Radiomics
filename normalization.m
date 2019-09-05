function norm_box = normalization(box)
max_value=max(box(:));
min_value=min(box(:));
for i=1:size(box,1)
    for j=1:size(box,2)
        for k=1:size(box,3)
            if box(i,j,k)==NaN
                norm_box(i,j,k)=NaN;
            else
                norm_box(i,j,k)=(box(i,j,k)-min_value)/(max_value-min_value);
            end
        end 
    end
end
end