
function job=starter(ii)

%ii=str2num(ii);

folder = 3020; % 30 variables 20 obs.

%Depending on folder
%number of endogenous variables
N = 30;
%time
T = 20;
%number of exogenous variables in each equation
n_x = 1;

% 1 Monte Carlo replcs each time
Mt = 1;

% coefficients, the last one is for spatial effects while the first two are
% for x_i of each individual equation
theta = [0.9; 0.6];

% number of spatially related units on each side
m = 2;

%% Priors

% 1st stage
n = N-1;   % number of variables in y but for y_now
p = n_x*N; % number of vairables in x_all
np = n*p;
nn2 = .5*(n^2-n);

%d_l for coefficients
a_11 =  1/np  ;
%covariance matrix
%exponetial for diagonal
s_10 = 0.01;
%d_l for off diagonal
a_12 = 1/nn2  ;

%2nd stage
n = N-1+n_x;
% d_l for coefficients
a_21 = 1/2 ;
% gamma for sigma
s_20 = 1;
niu_20 = 1;



job = main(folder, N, m, T, theta,n_x, Mt, a_11, a_12, s_10, a_21, s_20, niu_20,ii);
end