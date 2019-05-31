function [ W, real_mean, real_std ] = weight_maker(RanN,InpN,imean_weight,iwsig,dist_var )
%Create weights

switch dist_var
    case 'gaus' 
    W=imean_weight+iwsig.*randn(RanN,InpN);
    
    case 'powr'
        A=(2+sqrt(4+4*imean_weight^2/iwsig^2))/2;
        K=imean_weight*(A-1)/A;
        
        W=randraw('pareto', [K, A], [RanN, InpN]);
    
    case 'expo'
        W=exprnd(imean_weight,RanN, InpN);
        disp('STD cannot be controlled with an exp dist')
    case 'logn'
         MU = log(imean_weight^2 / sqrt(iwsig^2+imean_weight^2));
       SIGMA = sqrt(log(iwsig^2/imean_weight^2 + 1));
       W=lognrnd(MU, SIGMA, RanN,InpN);
       
    case 'unif'
        A=(imean_weight-.5*sqrt(12*iwsig^2));
        B=2*imean_weight-A;
        
        W = A + (B-A)*rand(RanN,InpN);
        

        
end

W(W<0)=0;
real_mean=mean(W(:));
real_std=std(W(:));

end

