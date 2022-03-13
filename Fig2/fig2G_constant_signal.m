function[P]=fig2G_constant_signal(bs,r,open,T,initial,on,BF,cof,co)


off=on/(1/BF-1); %%%Set the SDTF unbinding rate
cf1=cof;       
cf2=1/cf1;

%%% For each state, define the unwrapping and rewrapping rates considering
%%% cooperativity 
 for i=1:15
    DNAopen(i)=open*(cf1)^(i-1);
    DNAclose(i)=open*co*(cf2)^(i-1);
end



%DNAopen(end)=0;
DNAclose(end-1)=0;
DNAclose(end)=0;





%%%Define the transition probability matrix when the non-osillcatory signal (See SI Eq 2)

Qone=zeros(15,15);
Qtwo=zeros(15,15);
for i=1:14
    Qone(i,i+1)=DNAopen(i);     
    Qone(i+1,i)=DNAclose(i);
end


for i=1:14  
    Qtwo(i,i+1)=DNAopen(i);     
    Qtwo(i+1,i)=(1-exp(-(1/r)*(bs-i-1)^2))*DNAclose(i);
end





Q33=zeros(15,15);
Q4=zeros(15,15);
for i=1:15
    Q33(i,i)=on/2;
    Q4(i,i)=off;
end    


Qc=[Qone,Q33;Q4,Qtwo];

for i=1:30  
    Qc(i,i)=-sum(Qc(i,:)); 
end





mu=zeros(1,30);
mu(initial)=1;
P=mu*expm(Qc*T); %%% Compute the probability distribution at time T

end