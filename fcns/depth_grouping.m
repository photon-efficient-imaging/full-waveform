function [inds_new,vals_new] = depth_grouping(x)

inds = find(x);
flags = [1;diff(inds)];
%
inds_new = []; 
vals_new = [];
inds_group = [];
for k=1:length(inds)
    if(flags(k)==1)
        inds_group = [inds_group,inds(k)];
        if(k==length(inds))
            inds_new = [inds_new; mean(inds_group)];
            vals_new = [vals_new; mean(x(inds_group))];
        end
    elseif(flags(k)~=1)
        inds_new = [inds_new; mean(inds_group)];
        vals_new = [vals_new; mean(x(inds_group))];
        inds_group = inds(k);
    end
end

end

