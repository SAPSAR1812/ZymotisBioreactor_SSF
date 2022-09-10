%% Constantsu
global Cpa Cpb Cps H L f lambda rhob rhoa rhos Vz kb Xm Y e Ta h ka ks T0 Tw Topt X0 U 
Cpa=1180;
Cps=2500;
f=0.00246;
h=95;
H=0.35;
ka=0.0206;
ks=0.3;
L=0.03;
U=1;
Acc=1;
Ach=1;
Ta=35;
T0=35;
Tw=35;
Topt=38;
Vz=0.01;
X0=0.001;
Xm=0.125;
Y=8.366*10^6;
e=0.35;
lambda=2414300;
rhoa=1.14;
rhos=700;
rhob=e*rhoa+(1-e)*rhos;
kb=e*ka+(1-e)*ks;
Cpb=(e*rhoa*(Cpa+f*lambda)+(1-e)*rhos*Cps)/rhob;
a=alpha();
b=beta();
g=gamma();
relax=0;
%% FTCS
Vz=0.01;
Tw=35;
l=10; %x
m=10; %z
n=100000; %t
tmax=8000;
deltax=1/l;
deltaz=1/m;
deltat=tmax/n;
A=zeros([l,m,n]);
A(:,1,:)=1-T0/Ta;
A(2,:,:)=A(1,:,:);
X=zeros([l,m,n]);


for k=2:n-1
    for i=1:l
       X(i,1,k+1)=X(i,1,k)+deltat*(H/Vz)*mu(A(i,1,k))*(X(i,1,k)+X0/Xm)*(1-X(i,1,k)-X0/Xm);
    end
    for j=1:m
       X(1,j,k+1)=X(1,j,k)+deltat*(H/Vz)*mu(A(1,j,k))*(X(1,j,k)+X0/Xm)*(1-X(1,j,k)-X0/Xm);
    end     
    
    for i=2:l-1
        
        for j=2:m
            mu1=mu(A(i,j,k));
            X(i,j,k+1)=X(i,j,k)+deltat*(H/Vz)*mu1*(X(i,j,k)+(X0/Xm))*(1-X(i,j,k)-(X0/Xm));
            
            if(j==m)
                A(i,j,k+1)=A(i,j,k)-deltat*(a/deltaz)*(A(i,j,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1));
            else
                A(i,j,k+1)=A(i,j,k)-deltat*(a/(2*deltaz))*(A(i,j+1,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1));
            
            end
            
        end
    end
    A(l,:,k+1)=A(l-1,:,k+1)+0.01*deltax*phi()*(A(l-1,:,k+1)-(T0/Ta)+(Tw/Ta));
    %A(:,:,k+1)=relax*A(:,:,k+1)+(1-relax)*A(:,:,k);
end


T=A.*Ta+T0;
Xreal=X.*Xm+X0;
% x profile for z=0.5H, times 1s, 10s, 100s

%% Reversing Velocity

l=10; %x
m=10; %z
n=100000; %t
tmax=8000;
deltax=1/l;
deltaz=1/m;
deltat=tmax/n;
Ttest=zeros(1000,1);
A=zeros([l,m,n]);
A(:,1,:)=1-T0/Ta;
A(2,:,:)=A(1,:,:);
X=zeros([l,m,n]);
s=zeros(10,1);
for i=1:5000
    s(i)=100*i+1;
end
count=1;
for k=2:50000
   
    if any(s(:)==k)
        count=count*-1;
        
        A(:,:,k)=flip(A(:,:,k),2);
        %X(:,:,k)=flip(X(:,:,k),2);
        
    end
    A(:,:,k+1)=A(:,:,k);
    for i=1:l
       X(i,1,k+1)=X(i,1,k)+deltat*mu(A(i,1,k))*(H/Vz)*(X(i,1,k)+X0/Xm)*(1-X(i,1,k)-X0/Xm);
    end
    for j=1:m
       X(1,j,k+1)=X(1,j,k)+deltat*mu(A(1,j,k))*(H/Vz)*(X(1,j,k)+X0/Xm)*(1-X(1,j,k)-X0/Xm);
    end    
    for i=2:l-1
        for j=1:m
            if(count==1 && j==1)
                
            else
            mu1=mu(A(i,j,k));
            X(i,j,k+1)=X(i,j,k)+deltat*mu1*(H/Vz)*(X(i,j,k)+X0/Xm)*(1-X(i,j,k)-X0/Xm);
            if(j==1)
                A(i,j,k+1)=A(i,j,k)-deltat*(a/deltaz)*(A(i,j,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1));
            
            elseif(j==m)
                A(i,j,k+1)=A(i,j,k)-deltat*(a/deltaz)*(A(i,j,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1));
            else
                A(i,j,k+1)=A(i,j,k)-deltat*(a/(2*deltaz))*(A(i,j+1,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1));
            
            end
            A(i,j,k+1)=relax*A(i,j,k)+(1-relax)*A(i,j,k+1); 
            end
        end
    end
    A(l,:,k+1)=A(l-1,:,k+1)+0.01*deltax*phi()*(A(l-1,:,k)+(T0/Ta)-(Tw/Ta));
    if(count==1)
        Ttest(k)=A(5,10,k);
    else
        Ttest(k)=A(5,1,k);
    end
    Ttest(k)=max(max(A(:,:,k)));
end

T=A.*Ta+T0;
Xreal=X.*Xm+X0;
% x profile for z=0.5H, times 1s, 10s, 100s

%% Constant Cooling Techniques
Tw=25;
Ph=0.; %0.04
Pc=0.001; %0.001
l=10; %x
m=10; %z
n=100000; %t
tmax=8000;
deltax=1/l;
deltaz=1/m;
deltat=tmax/n;
A=zeros([l,m,n]);
A(:,1,:)=1-T0/Ta;
A(2,:,:)=A(1,:,:);
Ac=ones([m,n]);
Ac=((Ta-35)/35)*Ac;
X=zeros([l,m,n]);
Vz=0.01;
b=beta();
Vc=0.1;
for k=2:n-1
    
    
    for i=1:l
       X(i,1,k+1)=X(i,1,k)+deltat*mu(A(i,1,k))*(H/Vz)*(X(i,1,k)+X0/Xm)*(1-X(i,1,k)-X0/Xm);
    end
    for j=1:m
       X(1,j,k+1)=X(1,j,k)+deltat*mu(A(1,j,k))*(H/Vz)*(X(1,j,k)+X0/Xm)*(1-X(1,j,k)-X0/Xm);
    end     
    for i=2:l-1
        
        for j=2:m
            mu1=mu(A(i,j,k));
            X(i,j,k+1)=X(i,j,k)+deltat*mu1*(H/Vz)*(X(i,j,k)+(X0/Xm))*(1-X(i,j,k)-(X0/Xm));
            
            %FTCS with dissipation (second order term in x added) for Ac
            if(j==m)
                A(i,j,k+1)=A(i,j,k)-deltat*(a/deltaz)*(A(i,j,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1))-deltat*Ph*(mean(A(i,:,k))-Ac(i,k));
                %Ac(m,k+1)=Ac(m,k)-(deltat/deltax)*Vc*(Ac(m,k)-Ac(m-1,k))+Pc*deltat*(A(l,m,k)-Ac(m,k));
    
            else
                A(i,j,k+1)=A(i,j,k)-deltat*(a/(2*deltaz))*(A(i,j+1,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1))-deltat*Ph*(mean(A(i,:,k))-Ac(i,k));
                %Ac(j,k+1)=Ac(j,k)-(deltat/(2*deltax))*Vc*(Ac(j+1,k)-Ac(j-1,k))+(0.5*Vc*Vc*(deltat/deltax)^2)*(Ac(j+1,k)-2*Ac(j,k)+Ac(j-1,k))+Pc*deltat*(A(l,j,k)-Ac(j,k));
        
            end
             
            
            
        end
        
        
    end
    A(l,:,k+1)=A(l-1,:,k+1)+0.01*deltax*phi()*(A(l-1,:,k)-(T0/Ta)+(Tw/Ta));
    
    

    
end


T=A.*Ta+T0;
Xreal=X.*Xm+X0;
% x profile for z=0.5H, times 1s, 10s, 100s


%% Control Cooling
relax=0;
Tw=25;
Ph=0; %0.009
Pc=0.001; %0.08
Pc1=0;
l=10; %x
m=10; %z
n=100000; %t
tmax=8000;
deltax=1/l;
deltaz=1/m;
deltat=tmax/n;
A=zeros([l,m,n]);
A(:,1,:)=1-T0/Ta;
A(2,:,:)=A(1,:,:);
Ac=ones([l,n]);
%Ac=-0*Ac;
X=zeros([l,m,n]);
Vz=0.01;
b=beta();
Vc=0.1;
Kpp=300/32;
Kss=0.04;
%Taa=zeros(n,1);
%Taa(1)=32;
%Taa(2)=32;
Kdd=0.1;
lambdafilter=1;
Gc=tf([-Kpp,Kss],[1 10 200]);
%Gc=tf([0.1,0.1,0.1],[0.1 0.1 0.1 0.1 0]);
Twarray=zeros(n,1)+Tw;
for k=2:n-1
    
    for i=1:l
       X(i,1,k+1)=X(i,1,k)+deltat*mu(A(i,1,k))*abs(H/Vz)*(X(i,1,k)+X0/Xm)*(1-X(i,1,k)-X0/Xm);
    end
    for j=1:m
       X(1,j,k+1)=X(1,j,k)+deltat*mu(A(1,j,k))*abs(H/Vz)*(X(1,j,k)+X0/Xm)*(1-X(1,j,k)-X0/Xm);
    end
    
    for i=2:l-1
        
        for j=2:m
            mu1=mu(A(i,j,k));
            X(i,j,k+1)=X(i,j,k)+deltat*mu1*(H/Vz)*(X(i,j,k)+(X0/Xm))*(1-X(i,j,k)-(X0/Xm));
            
            
            if(j==m)
                A(i,j,k+1)=A(i,j,k)-deltat*(a/deltaz)*(A(i,j,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1))-deltat*Ph*(mean(A(i,:,k))-Ac(i,k));
            else
                A(i,j,k+1)=A(i,j,k)-deltat*(a/(2*deltaz))*(A(i,j+1,k)-A(i,j-1,k))+deltat*(b/(deltax^2))*(A(i+1,j,k)-2*A(i,j,k)+A(i-1,j,k))+(g)*(X(i,j,k)-X(i,j,k-1))-deltat*Ph*(mean(A(i,:,k))-Ac(i,k));
            
            end
            %A(i,j,k+1)=relax*A(i,j,k)+(1-relax)*A(i,j,k+1); 
        end
        %FTCS with dissipation (second order term in x added) 
        %Ac(i,k+1)=Ac(i,k)-(deltat/(2*deltax))*Vc*(Ac(i+1,k)-Ac(i-1,k))+(0.5*Vc*Vc*(deltat/deltax)^2)*(Ac(i+1,k)-2*Ac(i,k)+Ac(i-1,k))+Pc*deltat*(mean(A(i,:,k))-Ac(i,k));
     
    end
    %Ac(l,k+1)=Ac(l,k)-(deltat/deltax)*Vc*(Ac(l,k)-Ac(l-1,k))+Pc*deltat*(mean(A(l,:,k))-Ac(l,k));
    
     
    %[yout,tout]=lsim(Gc,reshape(A(10,10,1:k),k,1),linspace(0,1,k));
    %Taa(k+1)=35*yout(end)+35;
    
    %Ac(1,k+1)=-1;
    A(l,:,k+1)=A(l-1,:,k+1)+0.01*deltax*phi()*(A(l-1,:,k)-(T0/Ta)+(Tw/Ta));
    
    Tw=25-1*(A(5,10,k)*Ta);
    Twarray(k)=Tw;
end


T=A.*Ta+T0;
Xreal=X.*Xm+X0;


%%
figure;
plot(1:9,T(1:9,5,1),1:9,T(1:9,5,35),1:9,T(1:9,5,570),1:9,T(1:9,5,1430),1:9,T(1:9,5,2880));
legend('t=0s','t=1s','t=20s','t=50s','t=100s')
title('z=0.5H')
xlabel('x')
ylabel('T in Celsius')
% z profile for x=0.5L times 1s,10s,100s
figure;
plot(1:10,T(5,1:10,1),1:10,T(5,1:10,35),1:10,T(5,1:10,570),1:10,T(5,1:10,1430),1:10,T(5,1:10,2880));
title('x=0.5L')
legend('t=0s','t=1s','t=20s','t=50s','t=100s')
xlabel('z')
ylabel('T in Celsius')
% Time evolution for x=0.5L, h=0.5H
figure;
plot(1:n,reshape(T(5,5,:),n,1))
title('x=0.5L','h=0.5H')
xlabel('Time in 0.01s')
ylabel('Temperature')

% t=50s solid curve
figure;
surf(T(1:10,1:9,n/2))
xlabel('z')
ylabel('x')
zlabel('T')
title('t=50s')

% t=100s solid curve
figure;
surf(T(1:10,1:9,end))
xlabel('z')
ylabel('x')
zlabel('T')
title('t=ends')

% Biomass generation due to reaction
figure;
surf(Xreal(:,:,end/20))
xlabel('z')
ylabel('x')
zlabel('X')
title('t=10s, Vz=0.01')


%% Functions

function f=mu(theta)
global Ta T0
T=theta.*Ta+T0;
A=7.483*10^7;
B=1.3*10^47;
E1=70225;
E2=283356;
f=A*exp(-E1./(8.314*(T+273.16)))./(1+B*exp(-E2./(8.314*(T+273.16))));
end
function f1=alpha()
global Cpa f lambda Cpb rhob rhoa 
f1=rhoa*(Cpa+f*lambda)/(rhob*Cpb);
end
function f=beta()
global kb H Cpb rhob Vz L
f=kb*H/(Cpb*rhob*Vz*L*L);
end
function f=gamma()
global rhos Y e Xm Cpb rhob Ta
f=(rhos*Y*(1-e)*Xm)/(Cpb*rhob*Ta);
end
function f=phi()
global h L kb
f=(h*L)/kb;
end
function f=P()
global U 
f=U;
end
