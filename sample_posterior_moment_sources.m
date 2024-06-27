%NOTE: Run invert_strainrate_for_backsip.m first

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SECTION

%name of mat file to save samples in
filename = 'samples_Ms';

%number of MCMC samples
Nsamples = 10^4;

%source correlations will be ignored beyond the specified distance
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



GGGTGGG = GGG'*GGG;
Weight_terms = W'*W;
if variable_rake
    Weight_terms = Weight_terms + Gdamp_perp'*Gdamp_perp;  %note, weights built into W and Gdamp_perp and Gdamp
end
if include_moment
     Weight_terms = Weight_terms + Gdamp'*Gdamp;
end

Gsharp = (GGGTGGG + Weight_terms)\GGG'; 
Cov_s = Gsharp*Gsharp'; %covariance of slip deficit rate


pms = [pm_top;pm];
if variable_rake
    pms = [pms;pms];
end

stepsize = .5*sqrt(diag(Cov_s));
all_Ms = zeros(length(mhat),Nsamples);


%xy coordinates of all sources
source_xy = [pms(:,6:7);nodes;nodes;nodes];  %nodes repeated for three moment sources at each node
UB_all = [UB;inf*ones(size(nodes,1)*3,1)];  %infinite upper and lower bounds for moment sources
LB_all = [0*LB;-inf*ones(size(nodes,1)*3,1)];

for k=1:length(mhat)    
    
     %find sources that are stongly correlated with current source
     distance = sqrt( (source_xy(k,1)-source_xy(:,1)).^2 + (source_xy(k,2)-source_xy(:,2)).^2 );

    
     
    ind = ((abs(Cov_s(k,:))/Cov_s(k,k)) > 0.5) & distance'<correlation_distance;
    position = find(ind)==k;
    
  
    
    inv_Cov_small = inv(Cov_s(ind,ind));
    mhat_unbounded_small = mhat_unbounded(ind);
    m_small = mhat(ind);
    stepsize_small = stepsize(ind);
    UB_small = UB_all(ind);
    LB_small = LB_all(ind);

    Ms = zeros(length(m_small),Nsamples);
    m_prev = m_small;
    

    %now conduct MCMC sampling of m(ind) 
    cnt2=0;
    log_prob = -.5*(m_small-mhat_unbounded_small)'*inv_Cov_small*(m_small-mhat_unbounded_small);
    
    for j=1:Nsamples*length(m_small)
        
        
        cnt = mod(j,length(m_small))+1;
        m_small(cnt) = m_small(cnt)+2*(rand(1)-.5)*stepsize_small(cnt);
        
        
            
            
            
        if m_small(cnt)<LB_small(cnt) | m_small(cnt)>UB_small(cnt)
            
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
  

