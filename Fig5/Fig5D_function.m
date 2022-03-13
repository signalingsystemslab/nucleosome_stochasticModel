function[p,Q]=Fig5D_function(bs,T,DNAopen,DNAclose,on,BF,r)
co=(1/BF-1);

off=on/co;







initial=1;


%%%%%%%Q1 without SDTF
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



%%%%%%% We assumed that the binding rate (on) is decreasing toward zero
%%%%%%% when the DNA-Histone contact region is being close to dyad.
Q3=zeros(15,15);
Q4=zeros(15,15);
for i=1:15
    Q3(i,i)=on*(0.5+0.5/(exp(-1*(i-bs))+1));
    Q4(i,i)=off;
%     kon(i)=Q3(i,i);
%     BFvary(i)=off/(off+Q3(i,i));
end

Q=[Q1,Q3;Q4,Q2];
for i=1:30
    Q(i,i)=-sum(Q(i,:));
end
    







mu=zeros(1,30);
mu(initial)=1;
P=mu*expm(Q*T); %%%Prob T<tau

p=(P(15)+P(30));


end