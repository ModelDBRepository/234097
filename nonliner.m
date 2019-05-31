function [th_act] = nonliner(activity, thresh)
%thresh needs to be length(numcells), and then repped to be size of
%activity
acts=size(activity); %cells x conds x trials
threshes=repmat(thresh,[1,acts(2:3)]);

 
 th_act=1./(1+exp(-1.*(activity-threshes))); %log

th_act(th_act<0)=0;
end