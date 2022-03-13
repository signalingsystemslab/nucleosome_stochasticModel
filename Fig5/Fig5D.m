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




ON=zeros(M,1);
RateBF=zeros(M,1);
DS=10*BS-70*ones(1,M);
on=100;
co=(1/BF-1);
off=on/co;
for i=1:M
    bs=BSstate(i);
     ON(i)=on*(0.5+0.5/(exp(-1*(bs-5))+1)); %%% ON is kon varying in the distance between the SDTF binding site and the DNA opening status
    RateBF(i)=off/(ON(i)+off); %%% RateBF is the unbinding fraction
end




for i=1:M
    bs=BSstate(i);   
    [Prob1(i),Q]=Fig5D_function(bs,T,DNAopen,DNAclose,ON(i),RateBF(i),r);   
end


%%% symmetrize the eviction probability
p=0.5;


for i=1:M   
    Prob3(i)=Prob1(i)*p+Prob1(M+1-i)*(1-p);  %%% This is the eviction probability of the stochastic model  
end




  figure
plot(DS,Prob3,'b','linewidth',2);
xlabel('Relative distance to dyad(bps)')
ylabel('P(X(T)=14)')
xlim([-140,140])
ylim([min(Prob3)-0.1,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')




%%%% plot kon varying with respect to the distance between the SDTF binding site and the DNA opening status
on=100;
co=(1/BF-1);
off=on/co;
for i=1:15
    kon(i)=on*(0.5+0.5/(exp(-1*(i-5))+1)); 
end
    
figure
plot(0:1:14,kon,'linewidth',2)
ylim([40,100])
xlabel('Relative distance to dyad(bps)')
ylabel('Kon')
set(gca,'Fontsize',20,'fontname','Times New Roman')


