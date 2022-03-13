clear




%%%% Load the searched parameter
load('parameter_LPS.mat')
cof=parameter(1);
ratio=parameter(2);
BF=parameter(3);
r=parameter(4);
open=parameter(5);

T=240; %%% Total time duration



%%%% To plot the green curve, load the searched parameter with the TNF-NFkB
%%%% data
% load('parameter_nfkb_final.mat')
% cof=parameter(1);
% ratio=parameter(2);
% BF=0.3;  %%% We slightly change BF from the searched parameter with the TNF-NFkB
% r=parameter(4);
% open=parameter(5);




DNAopen=zeros(15,1);
DNAclose=zeros(15,1);
for i=1:15
    DNAopen(i)=open*cof^(i-1);
    DNAclose(i)=ratio*open*cof^(-i+1);
end
    DNAopen(end)=0;
    DNAclose(end-1)=0;
    DNAclose(end)=0;

BSstate=-12:1:28;
BS=-13:1:27;
M=numel(BS);





for i=1:M
    bs=BSstate(i);
    Prob1(i)=Fig6C_function(bs,T,DNAopen,DNAclose,BF,r,r);  
end



p=0.5;

for i=1:M   
    Prob3(i)=Prob1(i)*p+Prob1(M+1-i)*(1-p);
    Prob4(i)=Prob1(M+1-i);
end



DS=10*BS-70*ones(1,M);
load('Prob_LPS_sym.mat')
figure
plot(-100:20:100,Prob_LPS_sym,'r','linewidth',2,'displayname','Experiment')


hold on
plot(DS,Prob3,'b','linewidth',2,'displayname','Direct fitting with gradient descent');
xlabel('Relative distance to dyad(bps)')
ylabel('P(X(T)=14)')
xlim([-100,100])
ylim([min(Prob3)-0.1,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')
legend boxoff




