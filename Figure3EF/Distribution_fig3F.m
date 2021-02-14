clear
M=300; %%%Number of samples generated 
Acc=zeros(M,1);  %%% Accessibility under non-oscillatory signals
Acco=zeros(M,1); %%% Accessibility under oscillatory signals
bs=0;            %%% SDTF binding site
r=100000000;    %%%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
T=360;      %%% Time
on=100;   %%%Set the maximum SDTF binding rates
BF=0.3;  %%% Time fraction SDTF is unbound
co=5;    %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
h=1.1;  %%% the cooperativity constant 'h' 
 X=[0:1:14,0:1:14];
 
 
for i=1:M
    initial=floor(2*rand)+1;
    open=0.1+0.15*(rand-1/2);    %%% the initial unwrapping rate 'a_1'
    P=bs_prob_constant(bs,r,open,T,initial,on,BF,h,co);    
    Acc(i,1)=sum(X.*P);

    PP=bs_prob_oscillatory(bs,r,open,T,initial,on,BF,h,co);
    Acco(i,1)=sum(X.*PP);
end



figure
xbin=1:0.5:14;
histogram(Acc(:,1),xbin,'normalization','probability')
xlabel('Mean Accessibility')
ylabel('Probabilities')
 title(['On=',num2str(on/2),', constat SDTF'])
set(gca,'fontsize',20,'fontname','Times New Roman')     


figure
histogram(Acco(:,1),xbin,'normalization','probability')
xlabel('Mean Accessibility')
ylabel('Probabilities')
title(['On=',num2str(on),', oscillatory SDTF'])
set(gca,'fontsize',20,'fontname','Times New Roman')



%%% Coefficient of variations for each model
sqrt(var(Acc(:,1)))/mean(Acc(:,1))
sqrt(var(Acco(:,1)))/mean(Acco(:,1))
