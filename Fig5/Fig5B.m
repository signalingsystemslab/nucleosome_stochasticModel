clear

%%%% Load the searched parameters

load('parameter_nfkb_final.mat')
h=parameter(1); %%% cooperativity constant
ratio=parameter(2); %%% ratio between b_1 and a_1
BF=parameter(3); %%% The unbinding fraction
r=parameter(4); %%% range of the SDTF effect 
open=parameter(5); %%% a_1

T=240; %%% Total time duration

%%%% For plotting the non-cooperative case, set h=1


r1=1/30; %%% kon is set to be a gaussian curve as on*(1-affinity*exp(-((8-n)^2)*r1)) where n  is the SDTF binding site; 
affinity=0.5; %%% affinity will be multiplied to kon at the most buried site (n=7)



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

BForiginal=BF;
ON=zeros(M,1);
RateBF=zeros(M,1);
DS=10*BS-70*ones(1,M);
for i=1:M
    bs=BSstate(i);
     ON(i)=100;
    RateBF(i)=BForiginal;
end


for i=1:M
    bs=BSstate(i);   
    [Prob1(i),Q,openrate,bindingfraction]=Fig5B_function(bs,T,DNAopen,DNAclose,ON(i),RateBF(i),r,affinity,r1);    
    kon(i)=openrate;
    BFvary(i)=bindingfraction;   
end


%%% symmetrize the eviction probability
p=0.5;

for i=1:M   
    Prob3(i)=Prob1(i)*p+Prob1(M+1-i)*(1-p);  %%% This is the eviction probability of the stochastic model 
end



figure
plot(DS,Prob3,'b','linewidth',2);
xlabel('Relative distance to dyad(bps)')
ylabel('P(X(4hours)=14)')
xlim([-100,100])
ylim([min(Prob3)-0.1,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')


%%% plot of the binding rate kon against the SDTF binding site
figure
plot(0:1:14,kon(14:1:28),'linewidth',2)
xlim([0,14])
ylim([40,100])
xlabel('Relative distance to dyad(bps)')
ylabel('Kon')
set(gca,'fontsize',20,'fontname','Times New Roman')


