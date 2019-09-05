function [TS,TS_Si] = getTS(ROIonly,MaskOnly)
mat_in = ROIonly;
mask_in = MaskOnly;

mask = mask_in;

%mask= handles.mask_vol_for_TA{idx_img_set}{idx_voi};
%mask2 = mask(2:end-1, 2:end-1, 2:end-1);%previous
mask2 = mask;%new
% size = 1 ==> Surrounding 6 voxels
% (0,1,2) denotes (less than, equal to, greater than)
%texture_spectrum_matrix = zeros(size(mat_in)-[2 2 2]);%previous
texture_spectrum_matrix = zeros(size(mat_in));%new

for idx1 = 2:size(mat_in,1)-1
    for idx2 = 2:size(mat_in,2)-1
        for idx3 = 2:size(mat_in,3)-1
              texture_spectrum_matrix(idx1-1, idx2-1, idx3-1) = ...
                  (mat_in(idx1,idx2,idx3) == mat_in(idx1-1,  idx2,    idx3)  )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1-1,  idx2,    idx3)    )* 2 + ...
                 3 * ((mat_in(idx1,idx2,idx3) == mat_in(idx1+1, idx2,    idx3)  )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1+1, idx2,    idx3)    )* 2) + ...
                 9 * ((mat_in(idx1,idx2,idx3) == mat_in(idx1,     idx2-1, idx3)  )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1,     idx2-1, idx3)    )* 2) + ...
                27* ((mat_in(idx1,idx2,idx3) == mat_in(idx1,     idx2+1, idx3) )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1,     idx2+1,idx3)    )* 2) + ...
                81* ((mat_in(idx1,idx2,idx3) == mat_in(idx1,     idx2, idx3-1)  )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1,     idx2,    idx3-1) )* 2) + ...
               243*((mat_in(idx1,idx2,idx3) == mat_in(idx1,     idx2, idx3+1) )* 1 + (mat_in(idx1,idx2,idx3) > mat_in(idx1,     idx2,    idx3+1))* 2) + ...
               1; % +1 means zeros is also accounted for
              
        end
    end
end

texture_spectrum_matrix = texture_spectrum_matrix.*mask2;
                  
%texture_spectrum_Si = zeros(1, max(texture_spectrum_matrix(:)));
texture_spectrum_Si = zeros(1, 729);

for idx4 = 1:double(max(texture_spectrum_matrix(:)))
    texture_spectrum_Si(idx4) = length(find(texture_spectrum_matrix==idx4));
end

TS= texture_spectrum_matrix;
TS_Si = texture_spectrum_Si;
% 
% 
% 
% for idx_direction = 1:size(direction_mat,1)
%     max_lines = direction_mat(idx_direction, 2);
%     max_num_mat = direction_mat(idx_direction, 3);
%     run_length_matrix = zeros(max(mat_in(:)), max(size(mat_in)));%, 9);
%     for idx_slice = 1:max_num_mat
%         I = mat_in(:, :, idx_slice);
%         for idx_intensity = 1:max(mat_in(:))
%             mat_isintensity = (I==idx_intensity);
%             for idx_line = 1:max_lines
%                 vec_line = extract_vec(mat_isintensity, direction_mat(idx_direction, 1), idx_line);                
%                 now_pos = find_first_in_vec(vec_line);
%                 if now_pos>0
%                     while (now_pos<=length(vec_line))
%                         cut_vec = vec_line(now_pos:end);
%                         run_length_now = determine_run_length(cut_vec);
%                         run_length_matrix(idx_intensity, run_length_now) = run_length_matrix(idx_intensity, run_length_now)+1;
%                         cut_vec = vec_line(now_pos+run_length_now:end);
%                         first_idx = find_first_in_vec(cut_vec);
%                         if first_idx > 0
%                             now_pos = now_pos +run_length_now+first_idx-1
%                         else
%                             break;
%                         end
%                     end
%                 end
%             end
%             
%         end
%         %stop
%     end
%     run_length_matrix_all(:,:,idx_direction) = run_length_matrix;
%     if idx_direction == 4
%         mat_in = permute(mat_in, [3 2 1]);
%     end
%     if idx_direction == 7
%         mat_in = permute(mat_in, [3 2 1]);
%         mat_in = permute(mat_in, [1 3 2]);        
%     end
% end
% 
% run_length_matrix_all = sum(run_length_matrix_all, 3);
% 
% return;
% %
% % run_length_matrix = zeros(max(mat_in(:)), max(size(mat_in)));%, 9);
% %
% % % same slice, 90 degrees == towards the top
% % for idx_intensity = 1:max(mat_in(:))
% %     mat_isintensity = (I==idx_intensity);
% %     for idx2 = 1:size(mat_in,2)
% %         now_row = find_first_in_vec(mat_isintensity(:,idx2));
% %         if now_row>0
% %             while (now_row<=size(mat_in,1))
% %                 cut_vec = mat_isintensity(now_row:end, idx2);
% %                 run_length_now = determine_run_length(cut_vec);
% %                 run_length_matrix(idx_intensity, run_length_now) = run_length_matrix(idx_intensity, run_length_now)+1;
% %                 cut_vec = mat_isintensity(now_row+run_length_now:end, idx2);
% %                 first_idx = find_first_in_vec(cut_vec);
% %                 if first_idx > 0
% %                     now_row = now_row+run_length_now+first_idx-1;
% %                 else
% %                     break;
% %                 end
% %             end
% %         end
% %     end
% % end




% img_as_vec = img_in(:);
%
% index_mapping_mat = make_index_matrix(size(img_in));
