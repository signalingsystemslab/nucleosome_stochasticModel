function[p]=Fig6D_function(bs,T,DNAopen,DNAclose,BF,r1,r2)
co=(1/BF-1);
on=100;
off=on/co;


% if bs==8
%     on=0;
% end




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
    if i+1<=bs
    Q2(i+1,i)=(1-exp(-r1*(bs-i-1)^2))*DNAclose(i);
    elseif i+1>bs
    Q2(i+1,i)=(1-exp(-r2*(bs-i-1)^2))*DNAclose(i);
    end
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
mu(initial)=1;
P=mu*expm(Q*T); %%%Prob T<tau

p=(P(15)+P(30));


end