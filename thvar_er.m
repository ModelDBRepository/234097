function [ th_vars ] = thvar_er( th_var, W )

    cell_mean=squeeze(sum(W,2));
    th_vars=cell_mean.*th_var; 

end

