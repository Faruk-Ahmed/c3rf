function L = seg2label(seg, labels)

L = zeros(length(unique(labels)), 1);

unique_seg_labels = setdiff(unique(seg(:)),0);

for i = 1:length(unique_seg_labels)
    p_inds = (seg == unique_seg_labels(i));
    sp_inds = double(labels) .* p_inds;
    
    li = setdiff(unique(sp_inds(:)),0);
    for j = 1:length(li)
        L(li(j)) = unique_seg_labels(i);
    end
end


