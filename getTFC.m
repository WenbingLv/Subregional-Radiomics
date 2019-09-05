function TFC = getTFC(ROIonly,Maskonly)
mat_in = ROIonly;
mask = Maskonly;

%mask2 = mask(2:end-1, 2:end-1, 2:end-1);
mask2 = mask;
%out_mat = zeros(size(ROIonly)-[2 2 2]);
out_mat = zeros(size(ROIonly));
squeezed_matrix = zeros(3, 3, numel(out_mat));

counter = 0;
for idx1 = 2:size(mat_in,1)-1
    for idx2 = 2:size(mat_in,2)-1
        for idx3 = 2:size(mat_in,3)-1
            counter = counter+1;            
            squeezed_matrix(:,:,counter) = [ ...
                mat_in(idx1-1,idx2,idx3) mat_in(idx1,idx2,idx3) mat_in(idx1+1,idx2,idx3); ...
                mat_in(idx1,idx2-1,idx3) mat_in(idx1,idx2,idx3) mat_in(idx1,idx2+1,idx3); ...
                mat_in(idx1,idx2,idx3-1) mat_in(idx1,idx2,idx3) mat_in(idx1,idx2,idx3+1)];
        end
    end
end

code_in_vec = compute_texture_coding(squeezed_matrix);

counter = 0;
for idx1 = 1:size(out_mat,1)
    for idx2 = 1:size(out_mat,2)
        for idx3 = 1:size(out_mat,3)
            counter = counter+1;            
            out_mat(idx1, idx2, idx3) = code_in_vec(counter);            
        end
    end
end

out_mat = out_mat .* mask2;
TFC = out_mat;
%out_mat = reshape(code_in_vec, size(out_mat));
