function out_code = compute_texture_coding(vec_in_all)

% compute code matrix
mat_out_order =[
    1     1     1 %3
    1     1     2 %4
    1     1     3 %5
    1     2     2 %5
    1     1     4 %6
    1     2     3 %6
    2     2     2 %6
    1     2     4 %7
    1     3     3 %7
    2     2     3 %7
    1     3     4 %8
    2     2     4 %8
    2     3     3 %8
    1     4     4 %9
    2     3     4 %9
    3     3     3 %9
    2     4     4 %10
    3     3     4 %10
    3     4     4 %11
    4     4     4 %12
    ];
mat_out_vec = mat_out_order * [100 10 1]';

% generate codes
out_code = zeros(1,size(vec_in_all,3));
for idx = 1:size(vec_in_all,3)
    vec_in = vec_in_all(:,:,idx);
    for now_idx = 1:3
        a = vec_in(now_idx, 1);
        b = vec_in(now_idx, 2);
        c = vec_in(now_idx, 3);
        if b==c || b==a
            if c==a
                code(now_idx) = 1; %equal
            else
                code(now_idx) = 2; %flat -> hill
            end
        elseif a<b
            if b<c
                code(now_idx) = 3; % a<b<c
            else
                code(now_idx) = 4; % a<b>c
            end
        else
            if b>c
                code(now_idx) = 3; % a>b>c
            else
                code(now_idx) = 4; % a>b<c
            end
        end
    end
    code = sort(code) * [100 10 1]';
    out_code(idx) = find(mat_out_vec==code);
end

return;



