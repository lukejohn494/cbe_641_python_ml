%CBE641 Computational Project Final code
%start code
clc; clear all; format compact
T=2
aa=10000;
tic
 
k=20 
    x=k;
    y=k;
 
S=zeros(x,y); %randomly synthesize the spins of the structure
J=1;
h=1;
 
for m=1:1:x
    for n=1:1:y
bset=[-1,1];
pos=randi(length(bset));
S(m,n)=bset(pos);
end 
end 
%spins have been randomized
hm=zeros(x,y);
 
for i=1:1:aa*x^2 %Start the secondary loop; ends after a successful sweep
i;
    Sn=S;
v=[1:1:x];
s1=datasample(v,1); %randomly select a spin to change
s2=datasample(v,1);
Sn(s1,s2)=-1*Sn(s1,s2);
H=sum(sum(hm));
m=s1;
n=s2;
c=m-1;
d=m+1;
e=n-1;
f=n+1;
        if m-1==0
            c=x;
        end
        if m+1==x+1
            d=1;
        end
        if n-1==0
            e=y;
        end
        if n+1==y+1
            f=1;
        end
dE=-2*J*(S(m,n)*S(c,n)+S(m,n)*S(d,n)+S(m,n)*S(m,e)+S(m,n)*S(m,f));
if dE>=0
    W=1;
    S=Sn; %disp('Move was initially accepted')
elseif dE<0
    W=exp((dE)/T);
    r=rand;
    if r<W
        S=Sn;% disp('Move was accepted randomly from Boltzmann distribution')
    end   
end
 
ss=mean(mean(S));
E1=zeros(x,x-1);
E2=zeros(x-1,x);
for n=1:1:x-1
    if n==1
        E1(:,n)=(S(:,x)+S(:,2)).*S(:,1);
    elseif n>=2 & n<=x-1
    E1(:,n)=S(:,n).*S(:,n+1); 
    end 
end
for m=1:1:x-1
    if m==1
        E2(m,:)=(S(x,:)+S(2,:)).*S(1,:);
        elseif m>=2 & m<=x-1
    E2(m,:)=S(m, :).*S(m+1,:);
    end 
end
EE=mean(mean(E1))+mean(mean(E2));
sss(i)=ss;
eee(i)=ee;
 
end %end the loop to find new S
 
for k=1:1:aa
  ssss(k)=sss(k*x^2);
  eeee(k)=eee(k*x^2);
end
 
 
clc
disp('Congratulations, the Ising Ferromagnetic model has collected the data!!')
toc
k=[1:1:aa]
plot(k,ssss)
figure 2
plot(k,eeee)
