%Dirichlete VB, simultaneous equations Y=XB+E;



function [B, count] = HoS_SE2(Y, X, a, a_o, s0)

% dbstop if warning
% format long
 
[T,n] = size(Y);
p = size(X,2);
Tn = T*n;
np = n*p; %number of coefficients in vec(B)
n2n = (n^2-n)/2; %Total number of off-diag elements in Omega divided by 2




%%value to be reused
y = reshape(Y,Tn,1);
XX = X'*X;

np12 = (np+1)/2;
np12_0 = (n2n+1)/2;




%% Index associated with calculating Omega
ind_all = zeros(n-1,n);
for i = 1:n
       if i==1  
       ind = [2:n]'; 
      elseif i==n
       ind = [1:n-1]'; 
      else
       ind = [1:i-1,i+1:n]';
       end
       
       ind_all(:,i) = ind;
end

%%VB estimations
tic;
   
%% starting values

%B related

inv_tau = 100;
inv_lambda = ones(np,1)*100;


b_last = 100*ones(np,1);
bnV=10^(-4)*ones(np,1);

%Omega related

phi_sq=zeros(n);
phi2_sq=zeros(n);
phi_inv_sq=zeros(n);
      

phi_o=10*gamrnd(1,1,n2n,1); phi_o=phi_o/sum(phi_o); 
phi2_o=phi_o.^2; phi_inv_o=1./phi_o;  
    k0=0;
    for i = 1:n-1
              phi_sq(i,(1+i):n)=  phi_o((k0+1):(k0+n-i),1);
              phi2_sq(i,(1+i):n)=  phi2_o((k0+1):(k0+n-i),1);
              phi_inv_sq(i,(1+i):n)=  phi_inv_o((k0+1):(k0+n-i),1);
              k0=k0+n-i;
    end
%    phi_sq=phi_sq+phi_sq'+eye(n);
   niu_0=phi2_sq+phi2_sq'+eye(n); %% just for fun, no other reasons to use these as starting values
   lambda_0=phi_inv_sq+phi_inv_sq'+eye(n); 
   tau_0 = 10;


%more starting values
Sigma = diag(diag(cov(Y)))+1e-10*ones(n);
Omega = inv(Sigma);
bnV_o = 1e-10*ones(n);



check=10; % a number that is greater than 1 so the loop can start
count=0; %number of VB iterations

while check>1
    
    %% B related parameters and hyperparameters

   %q(inv_psi)
    inv_psi = 1/(1+inv_tau);

    %q(inv_niu)
    inv_niu =1./(1+inv_lambda);

    %q(inv_tau)
    inv_tau = np12/(inv_psi+.5*(sum(bnV.*inv_lambda)));
 
    %q(inv_lambda)
    inv_lambda = 1./(inv_niu+.5*inv_tau*bnV);
   

    %construct V_inv
    V_inv = diag(inv_tau*inv_lambda);
   
        
    %q(b), b is called 'beta' in the math notes
    V=inv(kron(Omega,XX)+V_inv);
    V_jj=diag(V);
    
    b=V*(kron(Omega,X')*y);
    bnV=(b.^2+V_jj);
    
    B = reshape(b,p,n);
    S_big = Y-X*B;
    S_big = S_big'*S_big;
    
    %calculate the impacts on variance covariance matrix related 
    TR_sq = zeros(n); 
    for i=1:n
        for j=i+1:n  
           V_now = V ((i-1)*p+1:i*p, (j-1)*p+1:j*p); 
           TR_sq (i,j)= trace(XX*V_now);     
        end
    end
    TR_sq = TR_sq + TR_sq';
    
    %% Omega related paramters and hyperparameters
    
    %q(tau_o)
    
    % omega_vector = Omega(tril(true(size(Omega)),-1));
    bnV_o_vec = bnV_o(tril(true(size(bnV_o)),-1));
 
    niu_vec  = niu_0(tril(true(size(niu_0)),-1));
    inv_niu_vec = 1./niu_vec;
    lambda_vec  = lambda_0(tril(true(size(lambda_0)),-1));
    inv_lambda_vec = 1./lambda_vec;

   %q(inv_psi_0)
    psi_0 = (1+1/tau_0);

    
    %q(inv_tau_0)
    tau_0 = (1/psi_0+.5*(sum(bnV_o_vec.*inv_lambda_vec)))/np12_0;
  
 
    %q(inv_lambda)
    lambda_00 = (inv_niu_vec+.5/tau_0*bnV_o_vec);
        

    %q(inv_niu_0)
    niu_00 =(1+1./lambda_00);

 
niu_0=zeros(n);
lambda_0=zeros(n);
 
H_inv_sq = zeros(n);
    
 k0=0;
    for i = 1:n-1
 
              niu_0(i,(1+i):n)=  niu_00((k0+1):(k0+n-i),1) ;
              lambda_0(i,(1+i):n)=  lambda_00((k0+1):(k0+n-i),1) ;

              
              H_inv_sq (i,(1+i):n)=tau_0*lambda_00((k0+1):(k0+n-i),1);
                       
              k0=k0+n-i;
    end
 
   niu_0=(niu_0+niu_0')/2 ; 
   lambda_0=(lambda_0+lambda_0')/2 ;     
   
   H_inv_sq = (H_inv_sq +  H_inv_sq')/2;
   
   %% draw Omega
  
   for i=1:n
      ind = ind_all(:,i);     
      Sigma_11 = Sigma(ind,ind); sigma_12 = Sigma(ind,i);
      sigma_22 = Sigma(i,i);
      
      h_12 = H_inv_sq(ind,i);
      TR_21 = TR_sq (ind,i);
      
     
      s_21 = S_big(ind,i); s_22 = S_big(i,i);
      
      V_now = V ((i-1)*p+1:i*p, (i-1)*p+1:i*p); 
      
      TR_diag = trace(XX*V_now);
      
      TR_sq(i,i)= TR_diag;
      
      gb = s_22+s0+ TR_diag;
      b1 = T/gb; % random gamma with shape=T/2, rate=gb/2
      
      inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
      
      C = inv(gb*inv_Omega_11+ diag(h_12)) ;
     
      b2 = -C*(s_21+TR_21);
      omega_12 = b2; omega_22 = b1 + b2'*inv_Omega_11*b2;   
      
     
      
      bnV_o_now  = omega_12.^2 + diag(C);
      bnV_o(i,ind) = bnV_o_now;
      bnV_o(ind,i) = bnV_o_now;
         
      %%  update Omega, Sigma
      Omega(i,ind) = omega_12; Omega(ind,i) = omega_12;
      Omega(i,i) = omega_22;
      temp = inv_Omega_11*b2;
      Sigma_11 = inv_Omega_11 + temp*temp'/b1;
      sigma_12 = -temp/b1; sigma_22 = 1/b1;
      Sigma(ind,ind) = Sigma_11; Sigma(i,i) = sigma_22;
      Sigma(i,ind) = sigma_12; Sigma(ind,i) = sigma_12;
      
      %
      
   end
%%

% % 
      
      Diff = abs(b-b_last);  %to calculate the |differences| between ELBOs
      b_last=b;
     
      if Diff<10^(-4)
          check=-10; 
          
      end
 count = count +1 ;     
   
end

        
        
    
    
  