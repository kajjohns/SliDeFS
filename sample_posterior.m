%NOTE: Run invert_strainrate_for_backsip.m first

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SECTION

%name of mat file to save samples in
filename = 'samples_Ms';

%number of MCMC samples
Nsamples = 10^4;

%slip rate correlations will be ignored beyond the specified distance
correlation_distance = 20;

%END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Conduct unbounded least squares inversion 
mhat_unbounded = GG\dd; 
   

%linear error propagation
if use_Gaussian_prior
    GGG = [G./repmat(Esig,1,size(G,2));Gprior];
else
    GGG = [G./repmat(Esig,1,size(G,2))];   
end

if variable_rake
    Gsharp = (GGG'*GGG+weight_smooth^2*[pm_smooth pm_smooth]'*[pm_smooth pm_smooth] + Gperp'*Gperp)\GGG'; %note, rake_weight is built into Gperp
else
    Gsharp = (GGG'*GGG+weight_smooth^2*pm_smooth'*pm_smooth)\GGG';
end
Cov_s = Gsharp*Gsharp'; %covariance of slip deficit rate


pms = [pm_top;pm];
if variable_rake
    pms = [pms;pms];
end

stepsize = .5*sqrt(diag(Cov_s));
all_Ms = zeros(length(mhat),Nsamples);



for k=1:length(mhat)    
    
     %find patches that are stongly correlated with current patch
     distance = sqrt( (pms(k,6)-pms(:,6)).^2 + (pms(k,7)-pms(:,7)).^2 );
    
     
    ind = ((abs(Cov_s(k,:))/Cov_s(k,k)) > 0.5) & distance'<correlation_distance;
    position = find(ind)==k;
    
  
    
    inv_Cov_small = inv(Cov_s(ind,ind));
    mhat_unbounded_small = mhat_unbounded(ind);
    m_small = mhat(ind);
    stepsize_small = stepsize(ind);
    UB_small = UB(ind);
    
    Ms = zeros(length(m_small),Nsamples);
    m_prev = m_small;
    

    %now conduct MCMC sampling of m(ind) 
    cnt2=0;
    log_prob = -.5*(m_small-mhat_unbounded_small)'*inv_Cov_small*(m_small-mhat_unbounded_small);
    
    for j=1:Nsamples*length(m_small)
        
        
        cnt = mod(j,length(m_small))+1;
        m_small(cnt) = m_small(cnt)+2*(rand(1)-.5)*stepsize_small(cnt);
        
        
            
            
            
        if m_small(cnt)<0 | m_small(cnt)>UB_small(cnt)
            
            accept = 0;
            
            
        else
            
            log_prob_trial = -.5*(m_small-mhat_unbounded_small)'*inv_Cov_small*(m_small-mhat_unbounded_small);
    
            %metropolis acceptance rule
      
            rat=exp(log_prob_trial-log_prob);
            if rat>1
               accept=1;
            else
               r=rand;
               if r<rat
                  accept=1;
               else
                  accept=0;
               end
            end
        end


        if accept==1
            log_prob = log_prob_trial;
            m_prev = m_small;
        else
            m_small = m_prev;
        end
        
         if cnt==1
             cnt2 = cnt2+1;
            Ms(:,cnt2) = m_small;
         end
         
    end
     
    
  
    
    
    all_Ms(k,:) = Ms(position,:);
    
    disp(['Completed ' num2str(k/length(mhat)*100) ' % of samples'])
end 


 save(filename,'all_Ms','-v7.3')
  

