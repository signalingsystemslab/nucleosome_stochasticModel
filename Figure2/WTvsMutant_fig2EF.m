clear
M=200; %%% Number of generated cell enviroments
K=50;  %%% Number of simulations for each cell




load('parameter_nfkb_final.mat') %%%Load the fitted paramters 

cofmean=parameter(1); %%% Mean of the cooperativity constant 'h'
comean=parameter(2);  %%% Mean of the ratio of the rewrapping rates to the unwrapping rates 'b_1/a_1'
BFmean=parameter(3);  %%% Mean of the time fraction SDTF is unbound 'BF'.
r=parameter(4);   %%% Mean of the SDTF effect range
open=parameter(5); %%% Mean of the initial unwrapping rate 'a_1'


X_osc=zeros(M,1); %%% Accessibility under the oscillatory signal
X_const=X_osc; %%% Accessibility under the non-oscillatory signal
X_fold=X_osc; %%% Fold change between X_osc and X_const

for i=1:M   %%%simulation for each cell

    Xo=zeros(K,1);
    Xc=zeros(K,1);
    Xf=Xo;
    bs=floor(8*rand); %%% randomize the SDTF binding site for each cell
   
    open=open+open*(rand-1/2); %%% randomize the unwrapping rate for each cell
   
    
    for k=1:K  %%%%%generate K samples
            
            %%% randomize the paramters
            range=1/r+2*(rand-1/2);   
            BF=BFmean+0.2*(rand-1/2);
            co=comean+1*(rand-1/2);
            cof=cofmean+0.2*(rand-1/2);
           
            
            %%% Rund the Markov mocel
            [p,p_oc,p_c]=RunMarkovChain(range,co,BF,open,bs,cof); 
            Xf(k)=p;    %%% indicates the fold change 
            Xo(k)=p_oc; %%% indicates the nucleosome accessibility under the oscillatory signal
            Xc(k)=p_c;  %%% indicates the nucleosome accessibility under the non-oscillatory signal
    end
    
    X_osc(i)=mean(Xo);
    X_const(i)=mean(Xc);
    X_fold(i)= mean(Xf);

end

mm=max(max(X_osc),max(X_const));
X_osc=X_osc/mm;  %%%normaize the accessibility
X_const=X_const/mm; %%%normaize the accessibility
           
%%%%Compute the coefficient of variations
 cv1=sqrt(var(X_osc))/mean(X_osc); 
 cv2=sqrt(var(X_const))/mean(X_const);
 xbin=0:0.05:1;


 %%% Plot the distribution of the accessibilities
figure
histogram(X_osc,xbin,'normalization','probability','displayname',['WT, cv=',num2str(cv1)])

hold on
histogram(X_const,xbin,'normalization','probability','displayname',['MM, cv=',num2str(cv2)])
legend
xlim([0.2,1])
set(gca,'Fontsize',20,'fontname','Times New Roman')

%%%plot the fold change distribution
figure
meanfold=mean(X_fold);
fbin=0.5:0.01:1;
histogram(X_fold,fbin,'normalization','probability','displayname',['Fold change, mean=',num2str(meanfold)])
set(gca,'Fontsize',20,'fontname','Times New Roman')

mean(X_fold)
