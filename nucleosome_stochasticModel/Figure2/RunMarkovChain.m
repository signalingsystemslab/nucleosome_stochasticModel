function[p,p_oc,p_c]=RunMarkovChain(range,co,BF,open,bs,cof)
%%%Set the maximum SDTF binding rates
on_osc=100;
on_const=100;
off=on_osc/(1/BF-1);  %%%Set the SDTF unbinding rate






%%% For each state, define the unwrapping and rewrapping rates considering
%%% cooperativity 
 for i=1:15
    DNAopen(i)=open*(cof)^(i-1);
    DNAclose(i)=co*open*(cof)^(-i+1);
end

DNAclose(end-1)=0;
DNAclose(end)=0;




%%%Define the transition probability matrix when the osillcatory signal is on the on-phase (See SI Eq 2)
Q1=zeros(15,15);
Q2=zeros(15,15);
for i=1:14
    Q1(i,i+1)=DNAopen(i);    
    Q1(i+1,i)=DNAclose(i);
end


for i=1:14
    Q2(i,i+1)=DNAopen(i);
    if i+1<=bs
    Q2(i+1,i)=(1-exp(-range*(bs-i-1)^2))*DNAclose(i);
    elseif i+1>bs
    Q2(i+1,i)=(1-exp(-range*(bs-i-1)^2))*DNAclose(i);
    end
end




Q3=zeros(15,15);
Q4=zeros(15,15);
for i=1:15
    Q3(i,i)=on_osc;
    Q4(i,i)=off;
end




Q=[Q1,Q3;Q4,Q2];
for i=1:30
    Q(i,i)=-sum(Q(i,:));
end


%%%initial probability 
mu=zeros(1,30);
mu(1)=1;


%%%Define the transition probability matrix when the osillcatory signal is on the off-phase (See SI Eq 2)
QQ=[Q1,0*Q3;Q4,Q2];
for i=1:30
    QQ(i,i)=-sum(QQ(i,:));
end

Q33=zeros(15,15);
for i=1:15
Q33(i,i)=on_const;
end


%%%Define the transition probability matrix under the non-osillcatory signal (See SI Eq 2)

Qc=[Q1,Q33;Q4,Q2];
 for i=1:30
   Qc(i,i)=-sum(Qc(i,:));
 end

    %%%From the data given in Figure 2A, we set the periods of the
    %%%oscillatory signal and the non-oscillatory signal are 50 and 240
    %%%minutes, respectively.
    %%%Hence the probability distribution can follows as shown in SI Eq 5
    %%%and 6.
    PP=mu*expm(Q*50)*expm(QQ*50)*expm(Q*50)*expm(QQ*50)*expm(Q*40);
    Pc=mu*expm(Qc*240);
    
    X=[0:1:14,0:1:14];
    p_oc=PP*X';
    p_c=Pc*X';
    
p=p_oc/p_c;

p_oc;
p_c;


end