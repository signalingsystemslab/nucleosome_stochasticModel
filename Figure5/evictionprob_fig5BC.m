clear


%%% Load previously searched parameter
 load('parameter_nfkb_final.mat')
%load('parameter_isre.mat')

h=parameter(1);   %%% the cooperativity constant 'h' 
ratio=parameter(2);  %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
BF=parameter(3);  %%% Time fraction SDTF is unbound
r=parameter(4);  %%%% SDTF effect range
open=parameter(5); %%% the initial unwrapping rate 'a_1'
T=240; %%% T=120 for the ISRE case


%%% For each state, define the unwrapping and rewrapping rates considering
%%% cooperativity 

DNAopen=zeros(15,1);
DNAclose=zeros(15,1);
for i=1:15
    DNAopen(i)=open*h^(i-1);
    DNAclose(i)=ratio*open*h^(-i+1);
end
    %DNAopen(end)=0;
    DNAclose(end-1)=0;
    DNAclose(end)=0;

    
    
    
BSstate=-12:1:28; %%% expanded state space
BS=-13:1:27;    
M=numel(BS);





for i=1:M
    bs=BSstate(i); %%% Binding site
    Prob1(i)=eviction(bs,T,DNAopen,DNAclose,BF,r); %%% Compute the full eviction probability   
end


%%% consider the symmetry (averaging out as describe in Fig 4E)
p=0.5;
for i=1:M   
    Prob2(i)=Prob1(i)*p+Prob1(M+1-i)*(1-p);    
end


DS=10*BS-70*ones(1,M);%% DS is the distance to the center from each state
plot(DS,Prob2,'b','linewidth',2);


xlabel('Relative distance to dyad(bps)')
ylabel('P(X(T)=14)')
xlim([-140,140])
ylim([min(Prob2)-0.1,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')




