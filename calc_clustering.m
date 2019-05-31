function [allPvals]=calc_clustering(activity_data, cue_mat)

%First do GLM analysis on each neuron to get coefs, then measure clustering

inputs=size(cue_mat,2);
act_stats=size(activity_data);

 if min(cue_mat(:,1))==1
 cue_mat(:,1)=cue_mat(:,1)-1;
 end




mod_warn=zeros(3,size(act_stats,1));
for c=1:act_stats(1) %for each cell
    time_mean=squeeze(activity_data(c,:,:));
    rp_nums=ones(inputs-1,1); rp_nums(1)=size(time_mean,2);

    
rep_cue_mat=repmat(cue_mat,rp_nums');
    cactiv=reshape((time_mean),act_stats(2)*rp_nums(1),1);
 mf1= GeneralizedLinearModel.fit(rep_cue_mat, cactiv, 'linear','Distribution', 'normal','Link','identity', 'CategoricalVars', [1:3],'DispersionFlag',true);
modlin_facts(c,:)=mf1.Coefficients.Estimate(2:end);
modlin_pvals(c,:)=mf1.Coefficients.pValue(2:end);

if size(lastwarn(),1)>0
    mod_warn(1,c)=1;
    lastwarn('')
end

end

sig_modlin_facts=modlin_facts;
sig_modlin_facts(modlin_pvals>.05)=0; %only use sig factors!
ust_struct=UniSphereTest(sig_modlin_facts,'test','bingham'); %THIS REQUIRES BRIAN LAU'S TEST OF UNIFORMITY ON THE HYPERSPHERE CODE. YOU CAN GET IT HERE: https://github.com/brian-lau/highdim 
allPvals(1)=ust_struct.pval; allPvals(2)=ust_struct.stat;

end

