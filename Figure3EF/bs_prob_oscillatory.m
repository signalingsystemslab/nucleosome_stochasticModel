function [P]=bs_prob_oscillatory(bs,r,open,Tmax,initial,on,BF,cof,co)

period=60;
N=floor(Tmax/period); %%Compute how many on-phase or off-phases should be there
Tr=Tmax-N*period;     %%Tr indicates the remaining time until Tmax after N on or off phases.
mu=zeros(1,30);    %%% mu is the initial probability 
mu(1)=1;
PP=eye(30);       






on=100; %%% Amplitude of the SDTF signal (SDTF binding rate)
off=on/(-1+1/BF);  %%% SDTF unbinding off rate


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
    Q2(i+1,i)=(1-exp(-(1/r)*(bs-i-1)^2))*DNAclose(i);
    
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

mu=zeros(1,30);
mu(1)=1;



%%%Define the transition probability matrix when the osillcatory signal is on the off-phase (See SI Eq 2)

QQ=[Q1,0*Q3;Q4,Q2];
for i=1:30
    QQ(i,i)=-sum(QQ(i,:));
end
    




if N==0  %%% if Tmax < Period
        PP=expm(Q*Tmax);
end
    
if N>0  %%% if Tmax > Period
  accumulated=0; %% indicates the time ellapsed
    for k=1:N 
        if Tmax-accumulated < period %%% if the remaining time < period
              break
        end
        
        if mod(k,2)==1        %%% the current time is on the on-phase
           PP=PP*expm(Q*period);
           accumulated=accumulated+period;
      
        else                 %%% the current time is on the off-phase
            PP=PP*expm(QQ*period);
            accumulated=accumulated+period;       
        end    
    end
   
    
    if Tr>0 && mod(N,2)==0     %%% N periods end and the remaining time is on the on-phase      
       PP=PP*expm(Q*Tr);     
    elseif Tr>0 && mod(N,2)==1 %%% N periods end and the remaining time is on the off-phase      
      PP=PP*expm(QQ*Tr);    
    end
    
end
%%Compute the probability distribution
P=mu*PP;
end