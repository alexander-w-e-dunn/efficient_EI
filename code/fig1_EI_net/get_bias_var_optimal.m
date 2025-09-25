
% computes bias of the estimators (E and I population) for the network
% optimizing the encoding error (type=1) and for the network optimizing the
% loss (type=2) or the error (type=1)

close all
clear
clc

% addpath([cd,'/code/function/'])
saveres=1;
showfig=1;

ntype={'min_error','min_loss'};
type=2;
disp(['computing bias of the estimator in network with ',ntype{type}]);

%% parameters

M=3;                                   % number of input variables    
N=400;                                 % number of E neurons   
nsec=2;                                % duration of the trial in seconds 

sigma_s=2;
tau_s=10;
tau_x=10;                              % time constant of the targets  

tau_e=10;                              % time constant of the excitatory estimates  
tau_i=10;                              % time const I estimates 

tau_re=10;                             % time const single neuron readout in E neurons
tau_ri=10;                             % time const single neuron readout in I neurons

% metabolic constant and noise strength for the networks optimizing error
% (type 1) and cost (type 2)
if type==1
    beta=6;
    sigmav=5.7;
else
    beta=14;
    sigmav=5;
end

dt=0.02;                               % time step in ms     
q=4;                                   % E-I ratio 
d=3;                                   % ratio of mean I-I to E-I connectivity

tau_vec=cat(1,tau_x,tau_e,tau_i,tau_re, tau_ri);

%% compute stimulus features and target signals

[s,x]=signal_fun(tau_s,sigma_s,tau_x,M,nsec,dt);

%% get decoding weights and recurrent connectivity

[w,J] = w_fun(M,N,q,d);

% AD: here I perturb the S matrix with different topology
% J is 4x1 cell, cell 1 is empty, 2,3 and 4 are II, IE, and EI connections

% get n.nodes
nnodes = size(J{2},1);
% get n. connections
nedges = size(J{2},1) * (size(J{2},2)-1);
% get mean weighted degree to preserve
round(mean(degrees_dir(J{2})));
% set density
dens = 0.9;
% set desired K
K = ceil(nedges*dens*1/size(J{2},1)); 
% J{2} = full(adjacency(WattsStrogatz(nnodes,nnodes,0.5)));
J{2} = impose_lattice_topology(J{2}, K);  % Convert J{2} (I-I connections) to lattice topology; use K nearest neighbors

%% compute performance 

T=nsec*1000/dt;
ntr=100;

estE=zeros(ntr,M,T);
estI=zeros(ntr,M,T);
tic
for ii=1:ntr
    
    [~,~,xhat_e,xhat_i] = net_fun_complete(dt,sigmav,beta,tau_vec,s,w,J);
    
    estE(ii,:,:)=xhat_e;
    estI(ii,:,:)=xhat_i;
    
end
toc

%% bias
% time-and dimension-dependent bias of the estimates

Bky=cell(2,1);
% AD: initialise Vky
Vky=cell(2,1);

bigx=permute(repmat(x,1,1,ntr),[3,1,2]);
Bky{1}=squeeze(mean(estE-bigx));
Bky{2}=squeeze(mean(estE-estI)); 
% AD: calculation of Vky missing so here:
Vky{1}=squeeze(var(estE-bigx));
Vky{2}=squeeze(var(estE-estI)); 

% prepare for boxplot
Bploty=cellfun(@(x) x(:),Bky,'un',0);
Vploty=cellfun(@(x) x(:),Vky,'un',0);

% average bias and variance
By=cellfun(@mean, Bploty);
Vy=cellfun(@mean, Vploty);

display(By,'average bias of estimators E and I')
display(Vy,'average variance of estimators E and I')

% prepare signals in 1 dimension for plotting
sigE=cell(2,1);
sigE{1}=squeeze(x(1,:));
sigE{2}=squeeze(mean(estE(:,1,:)));

sigI=cell(2,1);
sigI{1}=squeeze(mean(estE(:,1,:)));
sigI{2}=squeeze(mean(estI(:,1,:)));

% AD: tidx not computed so adding here in S
% tidx = (0:T-1) * dt/1000;   % time vector in seconds
tidx = (0:T-1) * dt;        % time in milliseconds

%%
if saveres==1
    
    param_name={{'N'},{'M'},{'tau_s'},{'beta'},{'sigmav'},{'tau_vec:X,E,I,rE,rI'},{'q'},{'d'},{'dt'},{'nsec'}};
    parameters={{N},{M},{tau_s},{beta},{sigmav},{tau_vec},{q},{d},{dt},{nsec}};
    
    % savefile='/Users/alexd/GitHub/efficient_EI/result/EI_net/';
    savefile='U:\GitHub\efficient_EI\result\EI_net\';
    savename=['bias_var_',ntype{type}];
    save([savefile,savename],'Bky','Bploty','Vky','Vploty','tidx','sigE','sigI','parameters','param_name')
end

%%

