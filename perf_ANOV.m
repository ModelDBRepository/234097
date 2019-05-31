function [selec_mats] = perf_ANOV(activity, cue_mat)
%creates a matrix of pairwise selectivity, the diagonal of which is the %
%of cells with each type of selectivity. One matrix for 2-way ANOVA and one
%matrix for 3-way

act_stats=size(activity);
inputs=size(cue_mat,2); 
for c=1:act_stats(1)
    time_mean=squeeze(activity(c,:,:)); 
    rp_nums=ones(inputs-1,1); rp_nums(1)=size(time_mean,2);
    rep_cue_mat=repmat(cue_mat,rp_nums');
    cactiv=reshape((time_mean),act_stats(2)*rp_nums(1),1);
    
  [pvs, tab, stats]=anovan(cactiv, rep_cue_mat, 'model', 'full','sstype', 2, 'display','off'); %3-way

  binp=double(pvs<.05); pures=(sum(binp(1:inputs))>0); mixs=(sum(binp(inputs+1:end))>0); binp(end+1)=pures; binp(end+1)=mixs; %total mix and pure
   bmat_full(c,:,:)=binp*binp';

  [pvs, tab, stats]=anovan(cactiv, rep_cue_mat, 'model', 'interaction','sstype', 2, 'display','off'); %2-way
  
 binp=double(pvs<.05); pures=(sum(binp(1:inputs))>0); mixs=(sum(binp(inputs+1:end))>0); binp(end+1)=pures; binp(end+1)=mixs;  
   
 bmat_inter(c,:,:)=binp*binp';

end




intermat=squeeze(sum(bmat_inter,1))./act_stats(1); fullmat=squeeze(sum(bmat_full,1))./act_stats(1); %percents
selec_mats.inter=intermat; selec_mats.full=fullmat;





end