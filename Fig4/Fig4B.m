clear

%%% Total time duration
 T=240;




%%%% Load the searched parameters

load('parameter_nfkb_final.mat')
h=parameter(1); %%% cooperativity constant
ratio=parameter(2); %%% ratio between b_1 and a_1
BF=parameter(3); %%% The unbinding fraction
r=parameter(4); %%% range of the SDTF effect 
open=parameter(5); %%% a_1



DNAopen=zeros(15,1);
DNAclose=zeros(15,1);
for i=1:15
    DNAopen(i)=open*h^(i-1);
    DNAclose(i)=ratio*open*h^(-i+1);
end
    DNAopen(end)=0;
    DNAclose(end-1)=0;
    DNAclose(end)=0;

BSstate=-12:1:28;
BS=-13:1:27;
M=numel(BS);





for i=1:M
    bs=BSstate(i);
    Prob1(i)=Fig4B_function(bs,T,DNAopen,DNAclose,BF,r,r);   
end


%%% symmetrize the eviction probability
p=0.5;

for i=1:M   
    Prob3(i)=Prob1(i)*p+Prob1(M+1-i)*(1-p); %%% This is the eviction probability of the stochastic model    
end


%%% Load the experimental data
DS=10*BS-70*ones(1,M);
load('Prob_nfkb_sym.mat')
figure
plot(-100:20:100,NFkB,'r','linewidth',2,'displayname','Experiment')


hold on
plot(DS,Prob3,'b','linewidth',2,'displayname','computation');
xlabel('Relative distance to dyad(bps)')
ylabel('P(X(T)=14)')
xlim([-100,100])
ylim([min(Prob3)-0.1,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')
legend boxoff



