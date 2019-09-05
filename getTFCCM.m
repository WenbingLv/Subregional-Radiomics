function [TFC,TFCCM] = getTFCCM(ROIonly,Maskonly)
mat_in = ROIonly;
mask = Maskonly;

%mask2 = mask(2:end-1, 2:end-1, 2:end-1);%previous
mask2 = mask;%new

%out_mat = zeros(size(ROIonly)-[2 2 2]);%previous
out_mat = zeros(size(ROIonly));%new
squeezed_matrix = zeros(3, 3, numel(out_mat));

counter = 0;
for idx1 = 2:size(mat_in,1)-1
    for idx2 = 2:size(mat_in,2)-1
        for idx3 = 2:size(mat_in,3)-1
            counter = counter+1;            
            squeezed_matrix(:,:,counter) = [
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
TFC  = out_mat;
%% co-occurrence matrix
img_in = out_mat; % Use the resampled image volume to compute the co-occurrence matrix
mask_in = mask2;

glm_3D = zeros(max(img_in(:))+1, max(img_in(:))+1); % intensity of zero is accounted as '1'

for idx1 = 1:size(img_in, 1)
    for idx2 = 1:size(img_in, 2)
        for idx3 = 1:size(img_in, 3)
            if mask_in(idx1, idx2, idx3)
                if idx2<size(img_in, 2) && mask_in(idx1, idx2+1, idx3)
                    a = sort([img_in(idx1, idx2, idx3) img_in(idx1, idx2+1, idx3)]); a = a+1;
                    glm_3D(a(1), a(2)) = glm_3D(a(1), a(2))+1;
                end
                if idx1<size(img_in, 1) && mask_in(idx1+1, idx2, idx3)
                    a = sort([img_in(idx1, idx2, idx3) img_in(idx1+1, idx2, idx3)]); a = a+1;
                    glm_3D(a(1), a(2)) = glm_3D(a(1), a(2))+1;
                end
                if idx3<size(img_in, 3) && mask_in(idx1, idx2, idx3+1)
                    a = sort([img_in(idx1, idx2, idx3) img_in(idx1, idx2, idx3+1)]); a = a+1;
                    glm_3D(a(1), a(2)) = glm_3D(a(1), a(2))+1;
                end
            end
        end
    end
end

for idx1 = 1:size(glm_3D,1)
    for idx2 = 1:size(glm_3D,2)
        if glm_3D(idx1, idx2) > 0
            glm_3D(idx2, idx1) = glm_3D(idx1, idx2);
        end
    end
end

% Create external global 'glcm_global' for other functions to access

glm_3D = glm_3D(2:end, 2:end);
TFCCM = glm_3D / sum(sum(glm_3D));

%out_mat = reshape(code_in_vec, size(out_mat));
