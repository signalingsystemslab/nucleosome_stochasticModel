clear
Tmax=500; %%%Maximum Time
bs=0 ; %%the NFkB binding site
range=100000; %%%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
count=0;



bs=0 ; %%the NFkB binding site
r=100000; %%%For Figure 3, we assume the SDTF effect range is from state 0 to state 14 (i.e. infinite effect range)
BF=0.3; %%% Time fraction SDTF is unbound
 co=15;   %%% the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
 cof=1.3;  %%% the cooperativity constant 'h'
 open=0.2; %%% the initial unwrapping rate 'a_1'
peri=0.005; %% Half-period of the SDTF signal



amplitude=100; %%% Amplitude of the SDTF signal (SDTF binding rate)
NFKBoff=amplitude/(-1+1/BF);    %%% SDTF unbinding off rate


%%% For each state, define the unwrapping and rewrapping rates considering
%%% cooperativity 
for i=1:15
    DNAopen(i)=open*(cof)^(i-1);
    DNAclose(i)=co*open*(cof)^(-i+1);
end
DNAopen(end)=0;
DNAclose(end-1)=0;
DNAclose(end)=0;





T=0:1:Tmax;

%%% Number of samples
K=50;



figure
parfor i=1:K
  
 
 
X=zeros(numel(T),1);       %%%Ttime evolution of DNA opening with m-th NFkB data.
x=0; %% state of DNA
t=0; %% Time 
nfkb=0; %% State of SDTF (either 0 or 1)
%%%%%%%%%%%%%%%%%%%%%%%%%% Gillespie simulation with WT NFkB dynamics
j=1;
while j<=numel(T)     

        
        
        
        
        
        
        %%%% Switch off backward jump at DNA=0. Absorb DNA at state 14
          
          
        
        
          if x==0
            switch1=0;
            switch2=1;
          elseif x==14
            switch1=0;
            switch2=0;
          else 
            switch1=1;
            switch2=1;
          end
        
        
        %%% Setting SDTF jump switch
          if nfkb==0
            switch3=0;
            switch4=1;
          else 
            switch3=1;
            switch4=0;
          end
          
          
          state=x;
          %%% SDTF effect
          if nfkb==1 
            switch5=0;
          else
            switch5=1;
          end
          
          
          
          
          
            n=0; 
                 while 1>0        %%% Oscillatory SDTF dynamics
                 if t< (2*n+1)*peri && t>=2*n*peri
                 r4=switch4*amplitude;
                 break
                 end
                 if t< (2*n+2)*peri && t>= (2*n+1)*peri
                 r4=0;
                 break
                 end
                 n=n+1;
                 end
        
        %%%%Setting Markov chain transition rates
        if x==0
            r1=0;
        else
        r1=switch1*DNAclose(x)*switch5;
        end
        
        r2=switch2*DNAopen(x+1);
        r3=switch3*NFKBoff;
        
               
        
        r0=r1+r2+r3+switch4*amplitude;
        
        tau=log(1/rand)/r0;
        t=t+tau;
       
        
         %%%%Store DNA when t > the corresponding time point  
     if t>=T(j)
         X(j,1)=x;
         j=j+1;
     end
       
      
     if x==14
         X(j:numel(T),1)=14;         
       break
     end
        
     q=j;
     %%% If t skips more than one time intervals, we fill X(j,m)=x for all
     %%% inbetween time points j
     for u=q:numel(T)
        if t >= T(u)
            X(u,1)=x;  
            j=j+1;
        else
            break
        end
     end
        
        
        
        
        
        
        
        
        
        
        r=rand;     
        
        if r<r1/r0
            x=x-1;      
        elseif r<(r1+r2)/r0
            x=x+1;        
        elseif r<(r1+r2+r3)/r0
            nfkb=0;
        elseif  r<(r1+r2+r3+r4)/r0
            nfkb=1;
        end
        
        
        
        
        
        
        if t>Tmax
             break
            end
            
     
     
     
     
     
     
end


%%% Check the trajectori reaches state 14
if x==14
    count=count+1;
    sw=1;
else 
    sw=0;
end


if sw==1
plot(T(1:numel(T)),X(1:numel(T)),'r','Linewidth',0.5)
else
plot(T(1:numel(T)),X(1:numel(T)),'k','Linewidth',0.5)
end
title(['DNA accessibility, period=',num2str(peri)])
ylabel('Chromatin accessibility')
xlim([0,Tmax])
ylim([0,15])
xlabel('Time (min)')
set(gca,'Fontsize',15,'fontname','Times New Roman')
hold on



end


