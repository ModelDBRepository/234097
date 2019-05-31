%% set parameters
clearvars
%----------------------network params
options=[2 4 4]; %number of task variable identity options per task variable 
remove_doubles=1; % if 1, don't allow 1st and 2nd cue to be the same

th_var=.27;  %threshold parameter

icon_prob=.25; %connection probability 


baseline=50; %baseline number of cells per population
incs_vec=[1.6,1,1.2]; %relative number of cells in each task variable variable population

  imean_weight=.207; %mean of weight distribution 
  iwsig_perc=1; %ratio of std of weight distribution to mean
  dist_var='gaus'; %distribution of weights. Options: 'gaus','powr','expo','logn','unif'
trials=10; %trials per condition

%-----------------learning params
 
  Nl=3; %Number of populations to increase weights to
  max_step=8; %total number of learning steps (1st step=random network)
    freelearn=1; %1 for free learning, 0 for constrained
  delr=.2; %learning rate


%-----------------Activity params
  mFR=4.9; %desired mean firing rate
  gnois_var=.88; %multiplicative noise parameter
  addsig_perc=5.8; %additive noise parameter

  Nuse=90; %number of "PFC" cells to use in analysis. total PFC population size is equal to total input population size
  


%------------necessary computations

inputs=length(options);
cue_perms=prod(options(1:end));
CellsOptsInps=ones(size(options))*baseline;
OptsInt=options;  
CellsOptsInps=ceil(CellsOptsInps.*incs_vec);



InpN=sum((OptsInt.*CellsOptsInps)); RanN=floor(InpN);  
rw_cells=[1:RanN]; %if I want to only have subset of cells do learning



 

% make matrix of input activity for all conditions
cue_mat=NaN(cue_perms, inputs);

for oi=1:length(options)

    n=options(oi);
    if oi<length(options)
        reps=prod(options(oi+1:end));
    else
        reps=1;
    end
 
    chunk_reps=cue_perms/(reps*n);
   
    for cri=1:chunk_reps
    for pn=1:n
        starti=1+(cri-1)*n*reps+(pn-1)*reps;
    cue_mat(starti:starti+reps-1,oi)=pn*ones(reps,1);
   
   
    end
    end
end
if remove_doubles
   
    lcue=cue_mat(:,end);
    samecue_inds=find(cue_mat(:,end-1)==lcue);
    cue_mat(samecue_inds,:)=[];
    cue_perms=size(cue_mat,1);
end
clear inp_starts; clear inp_ends
inp_starts(1)=1; ii=1; %find input boundaries

for j=1:sum(OptsInt)
    inp_ends(j)=inp_starts(j)+CellsOptsInps(ii)-1;
    inp_starts(j+1)=inp_ends(j)+1; 
    if j==sum(OptsInt(1:ii))
        ii=ii+1;
    end
end

addsig=imean_weight*addsig_perc;
iwsig=iwsig_perc*imean_weight; 

%% do learning and create activity for each condition

for rwi=1:max_step %through all learning steps

disp(['learning step: ',num2str(rwi)])
if rwi==1 %make initial random weight matrix


[W, imean_weight, iwsig]=weight_maker(RanN,InpN,imean_weight,iwsig,dist_var);
W=W.*double(rand(RanN,InpN)<=icon_prob); %connection probability 
W(W<0)=0; %only positive weights


threshes=thvar_er(th_var, W); %calculate thresholds for each cell
    
else %apply learning step
    
    [cx, cy]=find(W>0);
    if freelearn %increase input populations without regard to task variable
    for c1=1:length(rw_cells) 
       
        c=rw_cells(c1); 
        c_inps=cy(cx==c); 
        cinps_val=W(c,c_inps); 
        
        for ii=1:length(inp_starts)-1 %find max based on weights
            ind_vec=[((c_inps>=inp_starts(ii)).*(c_inps<=inp_ends(ii)))==1]; 
        hvec(ii)=sum(cinps_val(ind_vec));
        end
        [val(c,:), mi(c,:)]=sort(hvec,'descend'); %sort populations according to total input weight
        
        
        mi_inps=squeeze(mi(c,:));

        add_inds=[];
        for ci=1:Nl
            add_inds=[add_inds, [inp_starts(mi_inps(ci)):inp_ends(mi_inps(ci))]]; %create list of inputs to be increased
        end

       nadd_inds=1:InpN; nadd_inds=setdiff(nadd_inds, add_inds); 
        t_inp=sum(W(c,:));
        
        W(c,add_inds)=W(c,add_inds).*(1+delr); %increase weights
        
        ai_sum=sum(W(c,add_inds));
      
        W(c,:)=(W(c,:)./sum(W(c,:))).*(t_inp); %normalize
        
     
     
    end

    
   else %constrained learning

    for c=1:length(rw_cells) 
       
        c=rw_cells(c);
        c_inps=cy(cx==c); 


        cinps_val=W(c,c_inps); % allow rw to noise
           
        class1=1; clear classmax_ii classmax_val
        for oii=1:length(OptsInt) %go over all input types 
            class_end=class1+OptsInt(oii)-1;
            clear hvec
            ii1=0;
        for ii=class1:class_end %find max based on weights
            ii1=ii1+1;
            ind_vec=[((c_inps>=inp_starts(ii)).*(c_inps<=inp_ends(ii)))==1]; 
        hvec(ii1)=sum(cinps_val(ind_vec));
        
        end
        ii_inds=[class1:class_end]; 
        [val, mi]=sort(hvec,'descend');
        classmax_ii(oii)=ii_inds(mi(1)); 
        classmax_val(oii)=val(1);
        
        classmax_ii2(oii)=ii_inds(mi(2)); %if Nl>3, will need to know 2nd strongest input population from each task variable
        classmax_val2(oii)=val(2);

        class1=class_end+1; 
        
        
        end
        
        [valB, miB]=sort(classmax_val,'descend'); %order each input population with constraint that they are from different task variables
        [valB2, miB2]=sort(classmax_val2,'descend');
        mi_inps=[classmax_ii(miB), classmax_ii2(miB2)];
       
        
        add_inds=[];
        for ci=1:Nl
            add_inds=[add_inds, [inp_starts(mi_inps(ci)):inp_ends(mi_inps(ci))]];
        end
        
        
       nadd_inds=1:InpN; nadd_inds=setdiff(nadd_inds, add_inds); 
        t_inp=sum(W(c,:));
        
        W(c,add_inds)=W(c,add_inds).*(1+delr); %increase weights
        
        ai_sum=sum(W(c,add_inds));
      
        W(c,:)=(W(c,:)./sum(W(c,:))).*(t_inp); %normalize
             
     
    
end 

    end %end of learning type conditional 
end 


    W(W<0)=0; %theres no reason that any weights should've gone negative with learning, but just to be safe.
    
    

    W_all(rwi,:,:)=W;

                    

%% run simulation with current weights 

activity1=runRCN_track( W, threshes, cue_mat, addsig, RanN, InpN, OptsInt,trials,inp_starts, inp_ends); %unscaled activity


tnum=(mFR*.9)/mean(activity1(:)); %scale to rough area of desired mean, before multiplicative noise is added
activity2=activity1(:,:,:).*tnum; 


g_add=normrnd(activity2,activity2.*gnois_var); %add multiplicative noise
g_add(g_add<0)=0; %mult noise will likely lead to negative firing rates, so fix that

rinds=1:Nuse; 
addFR=mFR-mean(g_add(:)); %in case multiplicative noise didn't end up skewing things up enough

activity=g_add(rinds,:,:)+ max(addFR,0); 

 
selec_mats=perf_ANOV(activity,cue_mat);

sel_matfull(rwi,:,:)=selec_mats.full; sel_matinter(rwi,:,:)=selec_mats.inter; %using 3-way or 2-way ANOVA






 clustval=calc_clustering(activity, cue_mat);
 
 ClustVals(rwi)=clustval(2);


activity_all(rwi,:,:)=mean(activity,3);

meanz=mean(activity(:,:,:),3); varz=(std(activity(:,:,:),[],3)).^2;

for c=1:size(activity,1)
ffts(c)=nanmean(squeeze(varz(c,:))./squeeze(meanz(c,:)));
end

meanzz=mean(meanz,2); varzz=std(meanz,[],2).^2; RVs=(varzz./meanzz);


FRMs=mean(mean(activity,3),2);

FFAs_rwi(rwi,:)=[mean(RVs), mean(ffts), mean(FRMs)]; %RV, FF_T, and mean FR


end
%% plot properties over the course of learning

figure; % FF_T, RV, Mixed, Pure, Clust
subplot(2,3,1)
hold all
plot(FFAs_rwi(:,2))
plot([1,max_step],[2.8,2.8],'k:')
title('FF_T')
subplot(2,3,2)
hold all
plot(FFAs_rwi(:,1))
plot([1,max_step],[1.1,1.1],'k:')
title('RV')
subplot(2,3,3)
hold all
plot(sel_matfull(:,end,end))
plot([1,max_step],[.51,.51],'k:')
title('% Mixed')
subplot(2,3,4)
hold all
plot(sel_matfull(:,end-1,end-1))
plot([1,max_step],[.86,.86],'k:')
title('% Pure')
subplot(2,3,5)
hold all
plot(ClustVals)
plot([1,max_step],[186,186],'k:')
xlabel('Learning Steps')
title('Clustering Value')
subplot(2,3,6)
hold all
bar([0.6556, 0.3333, 0.5333, 0.2667, 0.2778, 0.1444, 0.0778],'k')
plot(diag(squeeze(sel_matfull(end,1:end-2,1:end-2))))
title('Selectivity at End') %full selectivity profile after the last learning step

