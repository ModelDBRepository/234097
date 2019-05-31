function [ R_activity ] = runRCN_track( W, threshes, cue_mat, addsig, RanN, InpN, OptsInt,trials,inp_starts, inp_ends) %OptsInt, cue_mat, CellsOptsInps, icon_prob, RI_perc, nonlin,th_var, W, imean_weight, iwsig, noise_on)
%pass in 

ttypes=size(cue_mat,1); %trial types
activity=NaN(RanN,ttypes, trials);
addnoise=addsig.*randn([RanN,ttypes,trials]);



for tti=1:ttypes
    inp_act=zeros(InpN,1);
    for opt_i=1:length(OptsInt)
        j=cue_mat(tti,opt_i)+sum(OptsInt(1:opt_i))-OptsInt(opt_i);
        inp_act(inp_starts(j):inp_ends(j))=1;
    end

    for tri_num=1:trials

        activity(:,tti,tri_num)=W*inp_act+squeeze(addnoise(:,tti,tri_num));

    end
end

R_activity=nonliner(activity, threshes);




end
