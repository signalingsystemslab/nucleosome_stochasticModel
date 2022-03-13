clear

bs=0 ; %%the NFkB binding site
r=100000; %%%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
BF=0.3; %%% Time fraction SDTF is unbound
 co=15;   %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
 cof=1.3;  %%% the cooperativity constant 'h'
 open=0.2; %%% the initial unwrapping rate 'a_1'
period=0.005; %% Half-period of the SDTF signal

T=0:1:500; %%% Time points
M=numel(T); 
Prob=zeros(M,1); %%% Prob stores the full DNA eviction probability for each time
parfor m=1:M
    t=T(m);
    [p]=fig2D_function(period,bs,r,open,t,BF,cof,co); %%%% prob_curve outputs the full eviction probability
    Prob(m)=p;   
end


figure
plot(T,Prob)
xlabel('Time (min)')
ylabel('Probability')
set(gca,'fontsize',20,'fontname','Times New Roman')

