clc;
clear;
%PARAMETERS%

data=load ('C:\MATLAB6p5\work\Huck\SSPA\CISDM96to04m.txt');   % load the data matrix y

data=data';
y=data(2:109,:);
id=data(1,:);

[n m]=size(y);       % n is the sample size of data matrix y
                     % m is the number of models of data matrix y

B=1000;              % number of bootstrapping
s_level=0.05;        % significant level
Q=0.9;               % the porobability of picking the following sample



%==========================================================================
%   calculate the maximum of the sample means
%==========================================================================
    y_mean=mean(y);                   %calculate the sample mean.  This will give you row(1*m) vector.
    max_y_mean=max(y_mean);           %the maximum of the sample means
                                      % or the non-standardized SPA statistis

    
%==========================================================================
%   SPA procedure starts from here
%==========================================================================

     if  max_y_mean<=0;     
           disp('Accepting all of the null hypptheses because all of the statistics are less than 0')
          %If all of the statistics is non-positive, the Step-SPA test just accepts the null. 
     else;                      
         disp('Testing procedure continues because one of the statistics is greater than 0')
         %If one of the statistics is positive, then the Step-SPA testing procedure continues.
         
   %================================================================
   % Calculate the covariance martix defined in Step-SPA test paper 
   %================================================================
         y_demean=y-ones(n,1)*y_mean;     % generate the de-meaned data
         y_var=y_demean'*y_demean/n;      % the variaince term    
            for i=1:n-1;
                sigma=(y_demean(1:(n-i),:))'*y_demean((i+1):n,:)+(y_demean((i+1):n,:))'*y_demean(1:(n-i),:);
                y_var=y_var+(  (( ( (n-i)/n )*(1-Q)^(i))   +  ((i/n)*(1-Q)^(n-i)) ) * sigma)/n;
            end
         % recall that the NW type estiamtor is 
         %\Omega_0+ \sum^{nw}_{i=1} ( 1-(i/(nw+1)))*(\Omega_1+ \Omega_1')
                 
 
   %==================================================
   % Calculate the standardized SPA statistic
   %==================================================
         std_vector=((diag(y_var)').^0.5);              %calculate the standard deviation vector of the models
         %std_vector=ones(1,m);                         
         %if you want a non-standardized version SPA test, you can set
         %std_vector=ones(1,m), instead of the std of the models ;
         sspa_statistics=n^(0.5)*(y_mean./std_vector);         % the vector of the standardized statistics of all models.
                 
                         
   %==================================================
   % Calculate the re-centering vector
   %==================================================
         mu=zeros(1,m);     % re-centering vector
         for i=1:m;
            if (n^(0.5)*y_mean(1,i)/std_vector(1,i))<=(-(2*log(log(n)))^(0.5));
                mu(1,i)= y_mean(1,i)/std_vector(1,i);             
                %in SPA test, if \bar{y_i}/sigma_i <= -(ln(ln(n)))^(0.5),
                %the recentering function for model i is \bar{y_i}/std_i %
                %otherwise, the recentering function=0
            end;
         end;
                  
   %==================================================
   % bootstrap procedure starts from here
   %==================================================              
         boot_means=zeros(B,m);  
         boot_sspa_statistic=zeros(B,1);           % the matrix for the bootstrapped statistics
         yy=[y_demean'  y_demean']';   
         % a 2n x m matrix of de-meaned y; we stack one on another 
              
         for b=1:B;
             ran_idx=floor(rand(n,1)*n)+1;    %the random index matrix
              
             pr=rand(n-1,1);                  %the probability matrix that will decide if we should 
                                              %get the next observation or do a random draw
              for j=2:n;                              
                  if pr(j-1,1) < Q;                
                      ran_idx(j,1)=ran_idx(j-1,1)+1;  
                      % if the value is less than Q, we take the next one for next period;                                
                      % that is, we re-define ran_idx(j,1)=ran_idx(j-1,1)+1. Then the 
                      % probability of picking the next index will be Q
                      % or we randomly pick one for next period.
                      % That is we do not change ran_idx(j,1).
                 end
             end
             x=yy(ran_idx,:);                         % X is the bth bootstrap (de-meaned) sample  
             x_mean=mean(x)./std_vector;              % X_mean is the mean of the bth bootstrap (de-meaned and standardized) sample
             recentered_x=x_mean+mu;                  % recentered_x is the recentered mean
             boot_means(b,:)=recentered_x;
         end
         
         
             %==================================================
             %Step-SPA procedure starts from here
             %==================================================
             
             k=1;                       %start with step 1
             
             reject_1=ones(1,m);        % if reject_1(1,j)=1, then model j is not rejected yet.  
                                        % if reject_1(1,j)=0, then model j is rejected.
                                        % the procedure start with every model is not rejected
                                        % Hence, reject_1=ones(1,m). 
                                         
             reject_2=ones(1,m);        % denote the rejected models after kth step
                                        
             step=0;                    % number of the steps in this test
             
             while k<=m;          %the maximum steps of this procedure is m
                   
                   for b=1:B;
                       boot_sspa_statistic(b,1)=n^(0.5)*max(boot_means(b,:).*reject_1);
                       % boot_means(b,:).*a1 
                       % for example, if model j is rejected before this step. then reject_1(j,1)=0.  Hence, the
                       % bootstraped statistics for model j will be 0.
                       % if model j is not rejected before this step. then reject_1(j,1)=1.  Hence, the
                       % bootstraped statistics for model j will be the bootstrapped mean (centered and standardized) of model j.
                   end
                   boot_sspa_statistic=sort(boot_sspa_statistic,1); 
                   sspa_critical= boot_sspa_statistic(floor((1-s_level)*B),1); 
                   %sspa_critical will be the critical value at this step
                 
                 for jj=1:m;
                     if   sspa_statistics(1,jj)>sspa_critical;
                          reject_2(1,jj)=0;
                          % if model jj is rejected at this step or at previous steps, then we
                          % set reject_2(1,jj)=0;  otherwise, reject_2(1,jj)=1.
                     end
                 end
                 step=step+1;   %the step count
                 if reject_2==reject_1; k=m+1;
                     % if reject_2==reject_1, that means, there is no model rejected at
                     % this stage, so the procedure has to stop.  Hence, we
                     % set k=m+1;
                 else k=k+1; reject_1=reject_2;
                     % otherwise, we go to next step and set k=k+1;
                     % and set reject_1=reject_2 for the next stage.
                     
                 end 
                 
         end
                 model_rejected=find(reject_2==0);
                 id_rejected=id(model_rejected);
                 %model_rejected'
                 id_rejected'
                 
                 step
                 % will report the models we reject
                 % if no model is rejected, will report  "Empty matrix: 0-by-1"
        
   end
 
