
clear


BF=0.3;  %%% Time fraction SDTF is unbound
co=5;     %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
h=1.2;  %%% the cooperativity constant 'h' 
open=0.05;  %%% the initial unwrapping rate 'a_1'

 for i=1:15
    DNAopen(i)=open*(h)^(i-1);
    DNAclose(i)=co*open*(h)^(-i+1);
end


on=100;  %%%% 
off=on/(1/BF-1);


%DNAopen(end)=0;
DNAclose(end-1)=0;
DNAclose(end)=0;





%%%Define the transition probability matrix (See SI Eq 2)
Qone=zeros(15,15);
Qtwo=zeros(15,15);
for i=1:14
    Qone(i,i+1)=DNAopen(i);     
    Qone(i+1,i)=DNAclose(i);
end


for i=1:14  
    Qtwo(i,i+1)=DNAopen(i);     
    Qtwo(i+1,i)=0*DNAclose(i); %%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
end

Q3=zeros(15,15);
Q33=Q3;
Q4=zeros(15,15);
for i=1:15
    Q3(i,i)=on;
    Q33(i,i)=on/2;
    Q4(i,i)=off;
end



mu=zeros(1,30); %%% Initial probability 
mu(1)=1;
M=101;     
Fold=zeros(M,1); %%% m-fold change to the unwrapping rates
maxMfold=30;  %%% Maximum fold change
Prob_os=zeros(M,1); %%%Full DNA eviction probability under the oscillatory signal
Prob_const=zeros(M,1); %%%Full DNA eviction probability under the non-oscillatory signal
 
for k=1:M
    m=0.1+(10-1)*(k-1)/(M-1); %%% m indicates the fold change
    Fold(k)=(m);
    Q1=m*Qone;
    Q2=m*Qtwo;  
    Q=[Q1,Q3;Q4,Q2]; %%%Define the transition probability matrix when the osillcatory signal is on the on-phase (See SI Eq 2)
    QQ=[Q1,0*Q3;Q4,Q2]; %%%Define the transition probability matrix when the osillcatory signal is on the off-phase (See SI Eq 2)
    Qc=[Q1,Q33;Q4,Q2]; %%%Define the transition probability matrix under the non-osillcatory signal (See SI Eq 2)

    for i=1:30
        Q(i,i)=-sum(Q(i,:));
        QQ(i,i)=-sum(QQ(i,:));
        Qc(i,i)=-sum(Qc(i,:));
    end

   
    
     
    PP=mu*expm(Q*60)*expm(QQ*60)*expm(Q*60)*expm(QQ*60)*expm(Q*60)*expm(QQ*60);
    Pc=mu*expm(Qc*360);
    Prob_os(k)=PP(15)+PP(30);
    Prob_const(k)=Pc(15)+Pc(30);   
end


figure
plot(Fold,Prob_os,'b','linewidth',2,'displayname','Oscillatory SDTF')
hold on
plot(Fold,Prob_const,'r','linewidth',2,'displayname','Constant SDTF')
hold on

set(gca,'fontsize',20,'fontname','Times New Roman')
legend('boxoff','Location','best')
xlabel('Fold change')
ylabel('Prob(X(T)=14)')

