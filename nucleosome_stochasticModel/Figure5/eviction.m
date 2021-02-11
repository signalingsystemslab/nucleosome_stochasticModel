function[p,Q]=eviction(bs,T,DNAopen,DNAclose,BF,r)
co=(1/BF-1); 
on=100; %%% Amplitude of the SDTF signal (SDTF binding rate)
off=on/co; %%% SDTF unbinding off rate





%%%Define the transition probability matrix when the non-osillcatory signal (See SI Eq 2)

Q1=zeros(15,15);
Q2=zeros(15,15);
for i=1:14
    Q1(i,i+1)=DNAopen(i);    
    Q1(i+1,i)=DNAclose(i);
end


for i=1:14
    Q2(i,i+1)=DNAopen(i);    
    Q2(i+1,i)=(1-exp(-r*(bs-i-1)^2))*DNAclose(i);   
end




Q3=zeros(15,15);
Q4=zeros(15,15);
for i=1:15
    Q3(i,i)=on;
    Q4(i,i)=off;
end

Q=[Q1,Q3;Q4,Q2];
for i=1:30
    Q(i,i)=-sum(Q(i,:));
end
    






%%% initial distribution
mu=zeros(1,30);
mu(1)=1;
P=mu*expm(Q*T); %%% Compute the probability distribution at time T

p=(P(15)+P(30)); %%% Full eviction probability


end