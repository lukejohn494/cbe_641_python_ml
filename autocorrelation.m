%CBE641 Computational Project Final code

%start code

clc; clear all; format compact

T=[2:0.01:2.6];

aa=1000;

tic

u=1
 
for k=[10 20 30 50] 

    x=k;

    y=k;

    s=zeros(length(T),aa*x^2);

    xv=zeros(length(T),aa*x^2);
 U=zeros(aa,x)    
GRR=ones(1,x)
for j=1:1:length(T) %randomize the matrix every time the temperature is changed


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

Sn=S;

v=[1:1:x];

s1=datasample(v,1); %randomly select a spin to change

s2=datasample(v,1);

Sn(s1,s2)=-1*Sn(s1,s2);

H=sum(sum(hm));

m=s1;

n=s2;

        if m-1==0

            c=x;

        else c=m-1;

        end

        if m+1==x+1

            d=1;

        else d=m+1;

        end

        if n-1==0

            e=y;

        else e=n-1;

        end

        if n+1==y+1

            f=1;

        else f=n+1;    

        end

dE=-2*J*(S(m,n)*S(c,n)+S(m,n)*S(d,n)+S(m,n)*S(m,e)+S(m,n)*S(m,f));

if dE>=0

    W=1;

    S=Sn;

    %disp('Move was initially accepted')

elseif dE<0

    W=exp((dE)/T(j));

    r=rand;

    if r<W

        S=Sn;   

       % disp('Move was accepted randomly from Boltzmann distribution')

    else 

        S=S;

       % disp('Move was not accepted')

    end   

end

if i==u*x^2;
    
U(u,:)=S(x/2,:); 

u=u+1
end



 


end %end the loop to find new S

for r=1:1:x

G(r)=mean(U(:,1).*U(:,r));

end

 

 

for r=1:1:x

GG(r)=G(r)-mean(mean([U(:,1),U(:,r)]))^2;

end

 

GRR(j,:)=GG;
u=1

end %Complete temperature loop
if k==30, GG30=GRR
else GG50=GRR
end


GRRR=transpose(abs(GRR(:,1:x/2)));

r=[1:1:x/2];

for u=1:1:length(T)

    fitt=fit(transpose(r),(GRRR(:,u)),'exp1');

    cory=coeffvalues(fitt);

    corlen(u)=1/-cory(2);

end

if k==10
    cormel1=corlen;
elseif k==20
    cormel2=corlen;
elseif k==30
    cormel3=corlen;
elseif k==50
    cormel4=corlen;
end



end%complete lattice size loop

clc

disp('Congratulations, the Ising Ferromagnetic model has collected the data!!')

toc

%export data to Excel
xlswrite('corlength.xlsx', [transpose(cormel1), transpose(cormel2), transpose(cormel3), transpose(cormel4)], 'A1');
