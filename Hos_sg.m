%Dirichlete VB, single equation, Horse shoe prior

%
function [b, count]=Hos_sg(y,X, a, s0, niu0)  
T=size(y,1);
np=size(X,2); %number of coefficients



%%value to be reused
T2niu02=T/2+niu0; %for sig^-2
XX=X'*X;
np12 = (np+1)/2;

%%VB estimations
tic;
   
%%starting values
inv_tau = 100;
inv_lambda = ones(np,1)*100;

sig2=var(y); s1=sig2*T2niu02;
b=ones(np,1)*10^(-15); 
b_last = 100*ones(np,1);
bnV=10^(-4)*ones(np,1);



check=10; % a number that is greater than 1 so the loop can start
count=0; %number of VB iterations
while check>1
    
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
        
    %q(b), b is called 'theta' in the math notes
   
    V=inv(XX/sig2+V_inv);
    V_jj=diag(V);

    
    b=V*(X'*y)/sig2;
    %b2=b.^2;
    bnV=(b.^2+V_jj);
   
    Diff=max(abs(b-b_last));
  
      if Diff<10^(-8)
          check=-10;  
      end
   b_last = b;

    %q(sig^(-2))
    %sum_bs=sum(bnV.*psii_inv.*phi2_inv*tau2_inv);
    s1=((y-X*b)'*(y-X*b)+trace(XX*V))/2+s0;
    sig2=s1/T2niu02;
    
    count=count+1;%Number of vb iterations
    
end

