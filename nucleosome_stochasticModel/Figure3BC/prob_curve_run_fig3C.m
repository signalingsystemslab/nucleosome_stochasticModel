clear

bs=0 ; %%the NFkB binding site
r=100000; %%%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
BF=0.3; %%% Time fraction SDTF is unbound
co=6;   %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
cof=1.1; %%% the cooperativity constant 'h'
open=0.1; %%% the initial unwrapping rate 'a_1'
period=180; %% Half-period of the SDTF signal

T=0:1:720; %%% Time points
M=numel(T); 
Prob=zeros(M,1); %%% Prob stores the full DNA eviction probability for each time
for m=1:M
    t=T(m);
    [p]=prob_curve(period,bs,r,open,t,BF,cof,co); %%%% prob_curve outputs the full eviction probability
    Prob(m)=p;   
end
figure
plot(T,Prob)
xlabel('Time (min)')
ylabel('Probability')
set(gca,'fontsize',20,'fontname','Times New Roman')

