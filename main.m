function job=main(folder, n_a,m, T_b, theta, n_x,Mt,  a_11, a_12, s_10, a_21, s_20, niu_20, i_hpc)

   %% Generate the data
    rng(i_hpc)
    % Generate W, the spatial matrix
    % W_n = Wgen1(n_a,m);
    load W_30.txt
    % Standardize W_n
    W_n = normw(W_30); 
    
    S_n=eye(n_a)-theta(2,1)*W_n; %spatial related, should be invertable;
    
  
    % start Monte Carlo
    for mt_s=1:Mt
       
         % generate x_nt, u_nt, and y_nt_hat, the reduced form version of y_nt;
                x_n=cell(1,n_a);
                u_n=cell(1,n_a);
                y_n_hat=zeros(T_b,n_a);
                
                for i=1:n_a                    
                    x_n_i=randn(T_b,n_x);
                    x_n_i_bt=x_n_i*theta(1,1);
                    x_n{1,i}=x_n_i;
                    u_n_i=0.1*randn(T_b,1);
                    u_n{1,i}=u_n_i;
                    y_n_hat(:,i)=x_n_i_bt+u_n_i;
                end
                
                % generate y_nt
                y_nt = zeros(T_b,n_a);
                
                for t=1:T_b
                    y_t = y_n_hat(t,:);
                    y_t = S_n\y_t'; %inv(S_n)*y_t';
                    y_nt(t,:) = y_t';
                end
                
                
                Lambda_sem = zeros(n_a-1,n_a); 
                b_sem = zeros(1,n_a);
%               sig2_sem=zeros(n_a,1);

                %estimate SEM using 2sls
                x_all = zeros(T_b,n_x*n_a); %the matrix of instruments
                for i=1:n_a
                    x_all(:,i)= x_n{1,i};
                end
 

tic;    
count10 = 0;
count20 = 0;
                for i=1:n_a
                    y_x=y_nt;
                    y_x(:,i)=[];
                        


                      %% Bayesian 2sls
                      % 1st stage
                      tic
                      [B, count1] = HoS_SE2(y_x, x_all,a_11, a_12, s_10);
                      Yhat = x_all*B;
                      toc

                      
                      % 2nd stage
                      y_now=y_nt(:,i);
                      x_now=[Yhat x_n{1,i}];
                      [b_now, count2]=Hos_sg(y_now,x_now,a_21, s_20, niu_20);    
                    
                    Lambda_sem(:,i)=b_now(1:(n_a-1),:); % the coeff of spatial effect;
                    b_sem(:,i)=b_now(n_a,:);      % the coeff of exogenous variables
                    
                    count10 = count10 + count1;
                    count20 = count20 + count2;
           
                 end
                
tm = toc;   % time used for 2sls            
save (sprintf('Lambda_%d_%d.out', folder,i_hpc ), 'Lambda_sem', '-ASCII') 
save (sprintf('b_%d_%d.out',  folder, i_hpc ), 'b_sem', '-ASCII')
save (sprintf('count1_%d_%d.out', folder,i_hpc ), 'count1', '-ASCII') %number of vb replic of 1st regression
save (sprintf('count2_%d_%d.out', folder,i_hpc ), 'count2', '-ASCII') %number of vb replic of 2nd regression
save (sprintf('tm_%d_%d.out', folder,i_hpc ), 'tm', '-ASCII')                 
     end  
            
        

        
     job = 100;   
 end
   
   
   
 


