clear
load('parameter_nfkb_final.mat') %%%Load the fitted paramters 
cof=parameter(1); %%% Mean of the cooperativity constant 'h'
co=parameter(2);  %%% Mean of the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
BF=parameter(3);  %%% Mean of the time fraction SDTF is unbound 'BF'.
r=parameter(4);   %%% Mean of the SDTF effect range
open=parameter(5); %%% Mean of the initial unwrapping rate 'a_1'


%%%Set the maximum SDTF binding rates
on_const=100;
off=on_const/(1/BF-1);  %%%Set the SDTF unbinding rate








%%% For each state, define the unwrapping and rewrapping rates considering
%%% cooperativity 
 for i=1:15
    DNAopen(i)=open*(cof)^(i-1);
    DNAclose(i)=co*open*(cof)^(-i+1);
end

DNAclose(end-1)=0;
DNAclose(end)=0;




 Time=0:0.1:240; %%%Duration of the NFkB singal
 Prob=zeros(numel(Time),1);  %%%Full eviction probability
 Aver_Prob=zeros(numel(Time),1);  %%%Average of the Full eviction probability with each binding site
    
for i=1:15
    bs=i; %%% Obtain the full eviction probability for each binding site i
 


%%%Define the transition probability matrix under the non-osillcatory signal (See SI Eq 2)
Q1=zeros(15,15);
Q2=zeros(15,15);
for i=1:14
    Q1(i,i+1)=DNAopen(i);    
    Q1(i+1,i)=DNAclose(i);
end


for i=1:14
    Q2(i,i+1)=DNAopen(i);
    Q2(i+1,i)=(1-exp(-((bs-i-1)^2)*r))*DNAclose(i);     
end




Q4=zeros(15,15);
for i=1:15  
    Q4(i,i)=off;
end



Q3=zeros(15,15);
for i=1:15
Q3(i,i)=on_const;
end
Qc=[Q1,Q3;Q4,Q2];
 for i=1:30
   Qc(i,i)=-sum(Qc(i,:));
 end

 
%%%initial probability 
mu=zeros(1,30);
mu(1)=1;
 
 
    
 for k=1:numel(Time)
    
    T=Time(k);
    %%%Probability distribution at time T   
    Pc=mu*expm(Qc*T);
    
    %%%% compute the full DNA eviction probability at time T       
    Prob(k)=Pc(15)+Pc(30);    
 end
 
 Aver_Prob=Aver_Prob+Prob/15;
 
end
 plot(Time,Aver_Prob,'b','linewidth',2)
 ylim([0,1])
 set(gca,'fontsize',20,'fontname','Times New Roman')
 ylabel('Prob of Full DNA eviction')
 xlabel('Duration of NFkB Signal (min)')
