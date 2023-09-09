
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The well-balanced algorithem preserves stationary shockwaves for two-layer
% shallow water equations.

%this is part of a paper titled "preserving stationary discontinuities in
%two-Layer shallow water equations with a novel well-balanced approach".

% majid akbari* 2023.
%* Department of Civil Engineering, University of Sistan and
%Baluchestan, Zahedan, Iran.
%Email: majidakbari@pgs.usb.ac.ir  ORCID: 0000-0002-8659-3753 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear;
for kr=1:4
for kn=1:6

nk=[50 100 200 400 800 3200]; % mesh
n = nk(kn)+3;
t=0;
a=((n-3)/20);
dx = 1/a;
g = 9.81;
CFL=0.5;% timestep

H2(1:n+2) = zeros+1+0;
U(1:n+2) = zeros+0.0; 
Z(1:n+2) = zeros+0.0;

Ul1(1:n+2) = zeros;
i=6*a:10*a;
     Z(i+1)=0.05.*((i/a-8+2)); %#ok<SAGROW>
i=10*a:14*a;    
    Z(i+1)=0.2-0.0125.*((i/a-12+2)).^2; %#ok<SAGROW>

% empty matrices
S1(1:n+2) = zeros;c1(1:n+2) = zeros;
S(1:n+2) = zeros;S3(1:n+2) = zeros;
T(1:n+2) = zeros;Z1(1:n+2) = zeros;Up(1:n+2) = zeros;
B(1:n+2) = zeros;M(1:n+2) = zeros;H(1:n+2) = zeros;Um(1:n+2) = zeros;
C(1:n+2) = zeros;N(1:n+2) = zeros;J(1:n+2) = zeros;Zp(1:n+2) = zeros;
D(1:n+2) = zeros;S(1:n+2) = zeros;L1(1:n+2) = zeros;L2(1:n+2) = zeros;
E(1:n+2) = zeros;u1(1:n+2) = zeros;u2(1:n+2) = zeros;S2(1:n+2) = zeros;Zm(1:n+2) = zeros;BB(1:n+2) = zeros;
k(1:n+2) = zeros;Fu1(1:n+2) = zeros;
h1(1:n+2) = zeros;h2(1:n+2) = zeros;u1(1:n+2) = zeros;u3(1:n+2) = zeros;

B1(1:n+2) = zeros;C1(1:n+2) = zeros;C2(1:n+2) = zeros;C3(1:n+2) = zeros;
abslambda2(1:n+2)=zeros;z(1:n+2)=zeros;h(1:n+2)=zeros;L(1:n+2) = zeros;S5(1:n+2) = zeros;JJ(1:n+2) = zeros;BB(1:n+2) = zeros;
A(1:n+2)=zeros;h3(1:n+2) = zeros;
P(1:n+2)=zeros;up(1:n+2)=zeros;um(1:n+2)=zeros;Q(1:n+2)=zeros;DT(1:n+2) = zeros;
C4(1:n+2) = zeros;F2(1:n+2) = zeros;
S8(1:n+2) = zeros;F4(1:n+2) = zeros;
U1(1:n+2) = zeros;zl1(1:n+2) = zeros;zl(1:n+2) = zeros;
SS(1:n+2) = zeros;
c2(1:n+2)=zeros;z1(1:n+2)=zeros;H1(1:n+2) = zeros;
Hp(1:n+2)=zeros;SR(1:n+2)=zeros;SRe(1:n+2)=zeros;
Hm(1:n+2)=zeros;SL(1:n+2)=zeros;SLe(1:n+2)=zeros;
Um(1:n+2)=zeros;Up(1:n+2)=zeros;rh(1:n+2)=zeros;rh2(1:n+2)=zeros;X(1:n+2)=zeros;T1(1:n+2)=zeros;Tl11(1:n+2)=zeros;
F(1:n+2) = zeros;
Fu(1:n+2) = zeros;

q1(1:n+2) = zeros;
Fl1p(1:n+2)=zeros;Fl1m(1:n+2) = zeros;q2(1:n+2) = zeros;

hl2(1:n+2)=zeros;hl3(1:n+2) = zeros;hhl2(1:n+2) = zeros;hhl3(1:n+2) = zeros;F3(1:n+2) = zeros;XX(1:n+2) = zeros;
SLe1(1:n+2)=zeros;SLe2(1:n+2)=zeros;SRe2(1:n+2)=zeros;SRe1(1:n+2)=zeros;

F4(1:n+2) = zeros;Sl1(1:n+2) = zeros;


Fu2(1:n+2)=zeros;
Fu2(1:n+2)=zeros;

FHpl1(1:n+2)=zeros;
FHml1(1:n+2)=zeros;
FUpl1(1:n+2)=zeros;
FUml1(1:n+2)=zeros;

Fup2(1:n+2)=zeros;
Fum2(1:n+2)=zeros;

Fp(1:n+2)=zeros;
Fm(1:n+2)=zeros;

FHp(1:n+2)=zeros;
FHm(1:n+2)=zeros;
FUp(1:n+2)=zeros;
FUm(1:n+2)=zeros;

% inital dam break
H(1:(n-2)/2)=0.6-Z(1:(n-2)/2);
H((n-2)/2:end)=0.21-Z((n-2)/2:end);% lower layer water height 
Hl1= zeros+0.7-H-Z;% upper layer water height 



rrn=[0.8 0.92 0.98 0.99]; % table of density ratio

rr=rrn(kr);
r2=1-rr;


Ffix=1;
tend=3000.0005;
if kr==1
 tend=200.0005; 
elseif kr==2   
tend=500.0005; 
elseif kr==3
tend=1000.0005;
elseif kr==4
tend=2000.0005;
end

while  t <= tend
    
     

   v1 = max(abs(Ul1./Hl1) +...
        sqrt(g*(Hl1)));
    v2= max(abs(U./H) +...
        sqrt(g*(rr.*Hl1+H)));

  v = max(v1,v2);
  dt = CFL*(dx/v); % timestep is chosen dynamically
 
   t = t+ dt;

% draw subplot 1

xx=(1:n-3)./a-dx;
if t>0
    
subplot(3,1,1)

x=(1:n-3)./a+0.0*dx;
i=2:n-2;
plot(x,(H(i)+Z(i)),'co-',x,(Hl1(i)+H(i)+Z(i)),'bo-',x,SS(i),'r^-','MarkerSize',7.5,'LineWidth',1.0)

hold on
plot(x,Z(i),'k-.','LineWidth',1.0)

xlabel('x [m]');
ylabel('h_1+h_2 +Z [m]');


hold off

grid off

subplot(3,1,2)

x=(1:n-3)./a+0.0*dx;
i=2:n-2;

plot(x,(U(i)),'gs-',x,BB(i),'r^-','MarkerSize',7.5,'LineWidth',1.0)
hold off

xlabel('x [m]');
ylabel('q_2 [m^2/s]');

subplot(3,1,3)
plot(x,(0.5.*1./g.*(U(i)./H(i)).^2+rr.*Hl1(i)+H(i)+Z(i)),'r>-','MarkerSize',7.5,'LineWidth',1.0)

xlabel('x [m]');
ylabel('E_2 [m]');

if t<=dt
hfig1=figure(1);
set(hfig1,'position',[5 5 700 1800])
end

end

top = title('');


drawnow
i = 2:n;
Um(i)=(U(i));
Up(i)=(U(i+1));

Hm(i)=(H(i));
Hp(i)=(H(i+1));


Zm(i)=(Z(i));
Zp(i)=(Z(i+1));

xj=(1:n+2)./a;


%%%
Fum2(i)=(U(i).^2./((H(i)))+0.5.*g.*((H(i)).^2)); % layer 2
Fup2(i)=(U(i+1).^2./((H(i+1)))+0.5.*g.*((H(i+1)).^2));

Fl1p(i)=(Ul1(i+1).^2./((Hl1(i+1)))+0.5.*g.*((Hl1(i+1)).^2)); % layer 1
Fl1m(i)=(Ul1(i).^2./((Hl1(i)))+0.5.*g.*((Hl1(i)).^2));

H1(i)=(((U(i).^2)./(2*(H2(i)-Z(i)).^2)+g.*(H2(i)))./g);
%%

J(i)=0.5.*(Ul1(i+1)+Ul1(i));

q01(i)=J(i).^2./(Hl1(i).*Hl1(i+1));
%%%%%%%%%%%%%%%%%%%%
J(i)=0.5.*(U(i+1)+U(i));

q02(i)=J(i).^2./(H(i).*H(i+1));
%%%%%%%%%%%%%%%%%%%%


B2(i)=0.5.*g.*(H(i)+H(i+1)).*r2.*(H(i+1)-H(i));

J(i)=0.5.*(U(i)+U(i+1));

Fu(i)=(J(i)).^2./H(i)+0.0.*(g.*H(i).^2);
Fu1(i)=(J(i)).^2./H(i+1)+0.0.*(g.*H(i+1).^2);

P(i)=-((Fu1(i)-Fu(i))+B2(i))./(H(i+1)-H(i));

%%%%%%%%%%%%
C0(i)=0.5.*(Hl1(i+1)+Hl1(i));
C01(i)=0.5.*(H(i+1)+H(i));

B02(i)=0.5.*g.*(H(i)+H(i+1)).*rr.*(Hl1(i+1)-Hl1(i));

B01(i)=0.5.*g.*(Hl1(i)+Hl1(i+1)).*(H(i+1)-H(i));

alpha=1;

Sl1s(i)=alpha.*(1.*-g.*(0.5.*Hl1(i+1)+0.5.*Hl1(i)).*(Z(i+1)-Z(i))+(q01(i)-g.*C0(i)).*(Hl1(i+1)-Hl1(i)));

S2s(i)=alpha.*(1.*-g.*(0.5.*H(i)+0.5.*H(i+1)).*(Z(i+1)-Z(i))-g.*0.5.*(H(i)+H(i+1)).*rr.*(Hl1(i+1)-Hl1(i)));

%%%%%%%%%%%%%%%%%%%%%% MLSWES  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% layer 1 %%%%%%%%%%%%%%%%%%%%%%%%%


Q1(i)=-(q01(i)-g.*C0(i)).*(Hl1(i+1)-Hl1(i));

F21(i)=(Q1(i))./(Hl1(i+1)-Hl1(i));

Q2(i)=-(q02(i)-g.*C01(i)).*(H(i+1)-H(i));

F42(i)=(Q2(i)-S2s(i))./(H(i+1)-H(i));

F22(i)=(+B01(i)-Sl1s(i))./(H(i+1)-H(i));

F41(i)=(+B02(i))./(Hl1(i+1)-Hl1(i));

P(i)=-(F21(i).*F42(i)-F22(i).*F41(i));
P(isinf(P))=0;
P(isnan(P))=0;


h(i)=H(i-1);
h1(i)=H(i-1);
h2(i)=H(i+1);

um(i)=(U(i-1))./H(i-1);
up(i)=(U(i+1))./H(i+1);

u1(i)=(Ul1(i-1))./Hl1(i-1);
u2(i)=(Ul1(i+1))./Hl1(i+1);

C1(i)=(P(i-1));
C2(i)=(P(i));
C3(i)=(P(i+1));


Q1(i)=-(q01(i)-g.*C0(i)).*(Hl1(i+1)-Hl1(i));

F21(i)=(Q1(i))./(Hl1(i+1)-Hl1(i));

Q2(i)=-(q02(i)-g.*C01(i)).*(H(i+1)-H(i));

F42(i)=(Q2(i))./(H(i+1)-H(i));

F22(i)=(+B01(i))./(H(i+1)-H(i));

F41(i)=(+B02(i))./(Hl1(i+1)-Hl1(i));

P(i)=-(F21(i).*F42(i)-F22(i).*F41(i));
P(isinf(P))=0;
P(isnan(P))=0;
C(i)=(P(i-1));
F(i)=(C(i-1));
C4(i)=(P(i+1));

%%% shock position
J(i)=min(C(i),max(P(i),C2(i)));
JJ(i)=min(F(i),max(C(i),C1(i)));

tt1=find(JJ>0 & J<0);

ttS=tt1;

%%
Fr1=H;
Fr2=H;
G=H;G1=H;
G2=H;G3=H;
c2(i)=(1-rr).*(g.*H(i)); %layer 2

c1(i)=(1-rr).*(Hl1(i).*g); %layer 1

Fr1(i)=(Ul1(i)./Hl1(i)).^2./c1(i); %layer 1

Fr2(i)=(U(i)./H(i)).^2./c2(i);%layer 2

G(i)=-(Fr1(i)+Fr2(i)-(1-rr).*Fr1(i).*Fr2(i)-1);
c2(i)=(g.*H(i)); %layer 2

c1(i)=(Hl1(i).*g); %layer 1

G(i)=((U(i)./H(i)).^2-c2(i)).*((Ul1(i)./Hl1(i)).^2-c1(i))-g.^2.*H(i).*Hl1(i).*rr;
G1(i)=G(i+1);
%%

tX1=find(B>0);


h(i)=Hl1(i+1);
ttY1=find(C>0 & P<0 & H>h2);
ttY2=find(P>0 & C4<0 & H<h2);


ttC1=find(C>0 & P<0);
ttC2=find(P>0 & C4<0);
ttC10=find(JJ>0 & J<0);
ttC20=find(C2>0 & C3<0);

ttE1=find(C>0 & P<0  & H<h2);
ttE2=find(P>0 & C4<0 & H<h2);

h(i)=H(i+2);

tYY1=find(C1>0 & J<0);
tYY2=find(J>0 & C3<0);
tYY10=find(C>0 & P<0);
tYY20=find(P>0 & C4<0);

ttX3=find(G>0 & G1<0 & h1>h2);


ttY=[ttY2 ttY1];
ttC=[ttC2,ttC1,ttC20,tt1];
ttE=[ttE1,ttE2];
tYY=[tYY1,tYY2,tYY10,tYY20];


th=1;find(h1<h2 & h>Hp);

ttD=ttY;

SS(:)=nan;
BB(:)=nan;


GG=H.*0+1;

GG(ttY1)=0;


if t>0
  %shock pos  
  disp("Xs=");disp(xx(ttS))  


SS(ttS)=H(ttS)+Z(ttS);
BB(ttS)=U(ttS);
end


FD=H;
FD(i)=1;
FD(ttC)=0;

FE=H;
FE(i)=1;
FE(ttE)=0;

FC(i)=1; %#ok<*SAGROW> 
FC(ttX3)=0;
%%

B2(i)=0.5.*g.*(H(i-1)+H(i+1)).*rr.*(Hl1(i+1)-Hl1(i-1));

c2(i)=sqrt(0.5.*g.*(H(i-1)+H(i+1)));
 S3(i)=-(c2(i).^2.*(Z(i+1)-Z(i-1)));
    
J(i)=0.5.*(U(i+1)+U(i-1));

Fu1(i)=(J(i)).^2./H(i+1)+0.5.*g.*(H(i+1)).^2;
Fu(i)=(J(i)).^2./H(i-1)+0.5.*g.*(H(i-1)).^2;



JJ(i)=-(Fu1(i)-Fu(i)-S3(i)+B2(i))./(H(i+1)-H(i-1));
JJ(isnan(JJ))=0;
JJ(isinf(JJ))=0;

%%%%%%%%%%%%%%%%%%%%%% MLSWES  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% layer 1 %%%%%%%%%%%%%%%%%%%%%%%%%

J(i)=0.5.*(Ul1(i+1)+Ul1(i-1));

qf1(i)=J(i).^2./(Hl1(i-1).*Hl1(i+1));
%%%%%%%%%%%%%%%%%%%%

J(i)=0.5.*(U(i+1)+U(i-1));

qf2(i)=J(i).^2./(H(i-1).*H(i+1));
%%%%%%%%%%%%%%%%%%%%%%

C(i)=0.5.*(Hl1(i+1)+Hl1(i-1));

C1(i)=0.5.*(H(i+1)+H(i-1));

Q1(i)=qf1(i)-(-g.*C(i)).*(Hl1(i+1)-Hl1(i-1));

Q10(i)=-(qf1(i)-g.*C(i)).*(Hl1(i+1)-Hl1(i-1));

Q2(i)=-(qf2(i)-g.*C1(i)).*(H(i+1)-H(i-1));

B2(i)=0.5.*g.*(H(i-1)+H(i+1)).*rr.*(Hl1(i+1)-Hl1(i-1));



Sl1s(i)=-g.*0.5.*(Hl1(i+1)+Hl1(i-1)).*(Z(i+1)-Z(i-1));

B1(i)=0.5.*g.*(Hl1(i+1)+Hl1(i-1)).*(H(i+1)-H(i-1))-Sl1s(i);

Sl1m(i)=-g.*0.5.*(Hl1(i-1)+Hl1(i-1)).*(Z(i+1)-Z(i-1));

B1m(i)=0.5.*g.*(Hl1(i-1)+Hl1(i-1)).*(H(i+1)-H(i-1))-Sl1m(i);

Sl1p(i)=-g.*0.5.*(Hl1(i+1)+Hl1(i+1)).*(Z(i+1)-Z(i-1));

B1p(i)=0.5.*g.*(Hl1(i+1)+Hl1(i+1)).*(H(i+1)-H(i-1))-Sl1p(i);

D(i)=(0.5.*(Ul1(i+1)+Ul1(i-1))-0.25.*(hl3(i+1)+hl2(i-1)).*(0.5.*(Ul1(i+1)+Ul1(i-1))./hl3(i+1)+0.5.*(Ul1(i+1)+Ul1(i-1))./hl2(i-1))).*(0.5.*(Ul1(i+1)+Ul1(i-1))./(hl3(i+1))-0.5.*(Ul1(i+1)+Ul1(i-1))./(hl2(i-1)));


S2s(i)=-g.*0.5.*(H(i+1)+H(i-1)).*(Z(i+1)-Z(i-1));


Sl1(i)=-g.*0.5.*(hl3(i+1)+hl2(i)).*(Z(i+1)-Z(i));
B1e(i)=(0.5.*g.*(hl3(i+1)+hl2(i)).*(H(i+1)-H(i)))-Sl1(i)-0.*D(i);

B1E(i)=B1(i)-D(i);((1-D(i)).*B1m(i)+D(i).*B1p(i));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F21s(i)=(Q10(i)-D(i))./(Hl1(i+1)-Hl1(i-1));

F42s(i)=(Q2(i)+B2(i)-S2s(i))./(H(i+1)-H(i-1));

F22s(i)=(Q10(i)+B1E(i))./(H(i+1)-H(i-1));

F41s(i)=B2(i)./(Hl1(i+1)-Hl1(i-1));


J(i)=-(F21s(i).*F42s(i)-F22s(i).*F41s(i));
J(isnan(J))=0;
J(isinf(J))=0;
%%%%%%%%%%%%%%%%%%%%%



X(:)=JJ(:);
J(:)=J(:);

X(J>0 & J>0)=0;
X(J<0 & J<0)=1;


J(i)=X(i);
JJ(i)=X(i);

z1(:)=Z(:);
z(:)=Z(:);
zl1=Z;
zl=Z;
h2=H;
h3=H;
hl2=Hl1;
hl3=Hl1;
hl02=Hl1;
hl03=Hl1;

h1(:)=H(:);
hl01(:)=Hl1(:);
Z1(:)=Z(:);


xj2(:)=xj(:);
xj3(:)=xj(:);
xj1(:)=xj(:);


if Ffix==0
tt1=[];
end

j=tt1;
js=ttS;

Z1(js)=J(js).*Z(js-1)+(1-J(js)).*(Z(js+1));

z1(js)=J(js).*Z(js)+(1-J(js)).*(Z(js+1));
z(js)=J(js).*Z(js-1)+(1-J(js)).*(Z(js));


h1(j)=J(j).*H(j-1)+(1-J(j)).*H(j+1);

xj1(j)=J(j).*xj(j-1)+(1-J(j)).*(xj(j+1));

h2(j)=J(j).*H(j)+(1-J(j)).*H(j+1);
h3(j)=J(j).*H(j-1)+(1-J(j)).*H(j);


u2(i)=U(i);
u3(i)=U(i);


j=tt1;


hl01(js)=JJ(js).*Hl1(js-1)+(1-JJ(js)).*Hl1(js+1);

hl02(js)=JJ(js).*Hl1(js)+(1-JJ(js)).*Hl1(js+1);
hl03(js)=JJ(js).*Hl1(js-1)+(1-JJ(js)).*Hl1(js);




%%

%%%%%%%%%%%%%%   MLSWEs fixer

Fu(i)=((U(i).^2)./h2(i))+0.5.*g.*(h2(i)).^2;
Fu1(i)=((U(i).^2)./h3(i))+0.5.*g.*(h3(i)).^2;

T2(i)=Fu1(i+1)-Fu(i);

Fu(i)=((U(i).^2)./h1(i))+0.5.*g.*(h1(i)).^2;
Fu1(i)=((U(i).^2)./h1(i))+0.5.*g.*(h1(i)).^2;

T1(i)=Fu1(i+1)-Fu(i);

Fu(i)=((U(i).^2)./h1(i))+0.5.*g.*((h1(i)).^2);


%%%% source term layer2

S2s(i)=-g.*0.5.*(h3(i+1)+h2(i)).*(z(i+1)-z1(i));


E(i)=Z(i+1)-Z(i);
E(E>0)=1;
E(E<0)=0;

S3(i)=(-g.*((1-E(i)).*(h3(i+1))+(E(i)).*(h2(i))-abs(z(i+1)-z1(i))./2).*(z(i+1)-z1(i)));

D(i)=(0.5.*(U(i+1)+U(i))-0.25.*(h3(i+1)+h2(i)).*(U(i+1)./h3(i+1)+U(i)./h2(i))).*(0.5.*(U(i+1)+U(i))./(h3(i+1))-0.5.*(U(i+1)+U(i))./(h2(i)));
D(isnan(D))=1;
D(isinf(D))=1;
D(D<0.0)=0;
D(D>1)=1;
%S2(i)=(1-D(i)).*S2s(i)+D(i).*(S3(i));


S2(ttY)=S2s(ttY);
S2(i)=-g.*0.5.*(h3(i+1)+h2(i)).*(z(i+1)-z1(i));
S2m(i)=-g.*0.5.*(h2(i)+h2(i)).*(z(i+1)-z1(i));
S2p(i)=-g.*0.5.*(h3(i+1)+h3(i+1)).*(z(i+1)-z1(i));

S1(i)=-g.*0.5.*(h1(i+1)+h1(i)).*(Z1(i+1)-Z1(i));


S1p(i)=-g.*0.5.*(h1(i+1)+h1(i+1)).*(Z1(i+1)-Z1(i));
S1m(i)=-g.*0.5.*(h1(i)+h1(i)).*(Z1(i+1)-Z1(i));
%%%%%%



Fu1(i)=(Ul1(i).^2./hl2(i))+0.5.*g.*((hl2(i)).^2);
Fu2(i)=(Ul1(i).^2./hl3(i))+0.5.*g.*((hl3(i)).^2);


Tl11(i)=Fu2(i+1)-Fu1(i);


%%%% sourcer term 2 layer 1
Sl1s(i)=-g.*0.5.*(hl3(i+1)+hl2(i)).*(zl(i+1)-zl1(i));

E(i)=Z(i+1)-Z(i);
E(E>0)=1;
E(E<0)=0;

S3(i)=(-g.*((1-E(i)).*(hl3(i+1))+(E(i)).*(hl2(i))-abs(zl(i+1)-zl1(i))./2).*(zl(i+1)-zl1(i)));

D(i)=(0.5.*(Ul1(i+1)+Ul1(i))-0.25.*(hl3(i+1)+hl2(i)).*(U(i+1)./h3(i+1)+Ul1(i)./hl2(i))).*(0.5.*(Ul1(i+1)+Ul1(i))./(hl3(i+1))-0.5.*(Ul1(i+1)+Ul1(i))./(hl2(i)))./(S3(i)-Sl1(i));


Sl1(ttY)=Sl1s(ttY);
Sl1(i)=-g.*0.5.*(hl3(i+1)+hl2(i)).*(zl(i+1)-zl1(i));
Sl1m(i)=-g.*0.5.*(hl3(i)+hl2(i)).*(zl(i+1)-zl1(i));
Sl1p(i)=-g.*0.5.*(hl3(i+1)+hl2(i+1)).*(zl(i+1)-zl1(i));
%%%%%



%%%%%%%%%%%%%%
D(i)=(0.5.*(U(i+1)+U(i))).^2.*(1./h1(i+1)-1./h1(i))-(0.5.*(h1(i+1)+h1(i))).*(0.5.*(U(i+1)+U(i))).^2.*(1./(2.*h1(i+1).^2)-1./(2.*h1(i).^2));



B21(i)=0.5.*g.*(h1(i+1)+h1(i)).*rr.*(hl01(i+1)-hl01(i))-0.*D(i);
B21e(i)=(0.5.*g.*(h1(i)+h1(i+1)).*rr.*(hl01(i+1)-hl01(i)))-S1(i)-0.*D(i);
B21p(i)=(0.5.*g.*(h1(i+1)+h1(i+1)).*rr.*(hl01(i+1)-hl01(i)))-S1p(i);


B21m(i)=(0.5.*g.*(h1(i)+h1(i)).*rr.*(hl01(i+1)-hl01(i)))-S1m(i);

D(i)=(-D(i)+B21e(i)-B21m(i))./(B21p(i)-B21m(i));
D(isnan(D))=0;
D(isinf(D))=0;
D(D<0.0)=0;
D(D>1)=1;

B21e(i)=FE(i).*((1-D(i)).*B21m(i)+D(i).*B21p(i))+(1-FE(i)).*B21e(i);

D(i)=(0.5.*(U(i+1)+U(i))).^2.*(1./h3(i+1)-1./h2(i))-(0.5.*(h3(i+1)+h2(i))).*(0.5.*(U(i+1)+U(i))).^2.*(1./(2.*h3(i+1).^2)-1./(2.*h2(i).^2));

D(ttD)=0;


B2(i)=0.5.*g.*(h3(i+1)+h2(i)).*rr.*(hl03(i+1)-hl02(i))+0.*D(i);
B2e(i)=(0.5.*g.*(h3(i+1)+h2(i)).*rr.*(hl03(i+1)-hl02(i)))-S2(i)-0.*D(i);
B2p(i)=(0.5.*g.*(h3(i+1)+h3(i+1)).*rr.*(hl03(i+1)-hl02(i)))-S2p(i);
B2m(i)=(0.5.*g.*(h2(i)+h2(i)).*rr.*(hl03(i+1)-hl02(i)))-S2m(i);

D(i)=(-D(i)+B2e(i)-B2m(i))./(B2p(i)-B2m(i));
D(isnan(D))=0;
D(isinf(D))=0;
D(D<0.0)=0;
D(D>1)=1;

B2e(i)=FE(i).*((1-D(i)).*B2m(i)+D(i).*B2p(i))+(1-FE(i)).*B2e(i);


D(i)=(0.5.*(Ul1(i+1)+Ul1(i))-0.25.*(hl3(i+1)+hl2(i)).*(0.5.*(Ul1(i)+Ul1(i+1))./hl3(i+1)+0.5.*(Ul1(i+1)+Ul1(i))./hl2(i))).*(0.5.*(Ul1(i+1)+Ul1(i))./(hl3(i+1))-0.5.*(Ul1(i+1)+Ul1(i))./(hl2(i)));


B1(i)=0.5.*g.*(hl3(i+1)+hl2(i)).*(H(i+1)-H(i))+0.*D(i);
Sl1(i)=-g.*0.5.*(hl3(i+1)+hl2(i)).*(zl(i+1)-zl1(i));
B1e(i)=(0.5.*g.*(hl3(i+1)+hl2(i)).*(H(i+1)-H(i)))-Sl1(i)-0.*D(i);
B1p(i)=(0.5.*g.*(hl3(i+1)+hl2(i+1)).*(H(i+1)-H(i)))-Sl1p(i);
B1m(i)=(0.5.*g.*(hl3(i)+hl2(i)).*(H(i+1)-H(i)))-Sl1m(i);

D(i)=(-D(i)+B1e(i)-B1m(i))./(B1p(i)-B1m(i));
D(isnan(D))=1;
D(isinf(D))=1;
D(D<0.0)=0;
D(D>1)=1;


B1e(i)=((1-D(i)).*B1m(i)+D(i).*B1p(i));

%%%%%%%%%%%%%%%%%%%%

U1(i)=(Ul1(i+1)+Ul1(i+1))./(Hl1(i)+Hl1(i+1));

U2(i)=(U(i-1)+U(i-1))./(H(i+1)+H(i));

V(i)=(U1(i)-U2(i)).^2./(g.*(1-rr).*0.5.*((Hl1(i)+Hl1(i+1))+(H(i)+H(i+1))));

%%%%%%%%%%%%%%%%%%%%
J(i)=0.5.*(Ul1(i+1)+Ul1(i));
q1(i)=J(i).^2./(hl2(i).*hl3(i+1));
%%%%%%%%%%%%%%%%%%%%
J(i)=0.5.*(U(i+1)+U(i));
q2(i)=J(i).^2./(h2(i).*h3(i+1));
%%%%%%%%%%%%%%%%%%%%
J(i)=0.5.*(U(i+1)+U(i));
q21(i)=J(i).^2./(h1(i).*h1(i+1));
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
J(i)=0.5.*(Ul1(i+1)+Ul1(i));
q01(i)=J(i).^2./(Hl1(i).*Hl1(i+1));
%%%%%%%%%%%%%%%%%%%%

J(i)=0.5.*(U(i+1)+U(i));

%J(tt3)=0.5.*(U(tt3+2)+U(tt3-1));

q02(i)=J(i).^2./(H(i).*H(i+1));
%%%%%%%%%%%%%%%%%%%%



B02(i)=0.5.*g.*(H(i)+H(i+1)).*rr.*(Hl1(i+1)-Hl1(i));

B01(i)=0.5.*g.*(Hl1(i)+Hl1(i+1)).*(H(i+1)-H(i));

C(i)=0.5.*(hl3(i+1)+hl2(i));

C1(i)=0.5.*(h3(i+1)+h2(i));

C21(i)=0.5.*(h1(i+1)+h1(i));


FD(i)=1;
FD(tYY)=0;

%%%%%%%%%%%%%%%%%%%%%% MLSWES  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% layer 1 %%%%%%%%%%%%%%%%%%%%%%%%%





Q01(i)=-(q01(i)-g.*C0(i)).*(Hl1(i+1)-Hl1(i));
Q02(i)=-(q02(i)-g.*C01(i)).*(H(i+1)-H(i));

Q2(i)=-(q21(i)-g.*C21(i)).*(h1(i+1)-h1(i));

F21(i)=Q01(i);

F42(i)=Q02(i);

F22(i)=+B01(i);

F41(i)=+B21(i);

Q1(i)=-(q1(i)-g.*C(i)).*(hl3(i+1)-hl2(i));
Q2(i)=-(q02(i)-g.*C01(i)).*(h1(i+1)-h1(i));


Q01(i)=-(q1(i)-g.*C(i)).*(H(i+1)-H(i)+Z(i+1)-Z(i));

F21(i)=FD(i).*Q1(i)+(1-FD(i)).*Q1(i);

F42(i)=Q02(i);

F22(i)=B01(i);

F21z(i)=FD(i).*Q1(i)+(1-FD(i)).*Q1(i);

F22z(i)=FD(i).*B01(i)+(1-FD(i)).*B01(i);

F41(i)=B02(i);

Q2(i)=-(q21(i)-g.*C21(i)).*(h1(i+1)-h1(i));

B210(i)=-0.5.*g.*(h1(i+1)+h1(i)).*rr.*(hl01(i+1)-hl01(i));

B210(i)=-0.5.*g.*(h1(i+1)+h1(i)).*rr.*(h1(i+1)-h1(i)+Z1(i+1)-Z1(i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q1(i)=-(q1(i)-g.*C(i)).*(hl3(i+1)-hl2(i));
Q2(i)=-(q2(i)-g.*C1(i)).*(h3(i+1)-h2(i));

F42s=H;
F21s=H;

F21s(i)=Q1(i);

F42s(i)=Q2(i)+B2e(i);

F22s(i)=Q1(i)+B1e(i);

F41s(i)=B2(i);
FS=F;
FS(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);

FS2(i)=F21(i).*F42(i)-F22(i).*F41(i);

F21s(i)=Q1(i);

F42s(i)=Q2(i)+B2e(i);

F22s(i)=Q1(i)+B1e(i);

F41s(i)=B2(i);

FS3=F;
FS3(i)=FS(i)./FS2(i);
FS3(isinf(FS3))=0;
FS3(isnan(FS3))=0;

%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q1(i)=-(q1(i)-g.*C(i)).*(hl3(i+1)-hl2(i));
Q2(i)=-(q21(i)-g.*C21(i)).*(h1(i+1)-h1(i));



F21s(i)=(Q1(i)+0.*B1(i)-0.*Sl1(i));

F42s(i)=(Q2(i)+B21e(i));

F22s(i)=(Q1(i)+B1e(i));

F41s(i)=(0.*Q2(i)+B21(i)+0.*S1(i));


FS(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);
FS2(i)=F21(i).*F42(i)-F22(i).*F41(i);
FS31=F;
FS31(i)=FS(i)./FS2(i);
FS31(isinf(FS31))=0;
FS31(isnan(FS31))=0;
%%%%%%%%%%%%%%%%%%%%%

FCE=H;
FCE(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);
%%%%%%%%%%%%%%%%%%%%%% MLSWEs
%%%%%%%%%%%%%%%%%%%%%% layer 2

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
Q1(i)=-(q1(i)-g.*C(i)).*(hl3(i+1)-hl2(i));
Q2(i)=-(q2(i)-g.*C1(i)).*(h3(i+1)-h2(i));


F21s(i)=(Q1(i)+B1e(i));

F42s(i)=(Q2(i)+0.*B2(i)-0.*S2(i));

F22s(i)=(0.*Q1(i)+B1(i)-0.*Sl1(i));

F41s(i)=(Q2(i)+B2e(i));

FS(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);
FS2(i)=F21(i).*F42(i)-F22(i).*F41(i);

F21s(i)=(Q1(i)+B1e(i));

F42s(i)=(Q2(i));

F22s(i)=(+B1(i));

F41s(i)=(Q2(i)+B2e(i));


FS4=F;
FS4(i)=FS(i)./FS2(i);
FS4(isinf(FS4))=0;
FS4(isnan(FS4))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% MLSWEs
%%%%%%%%%%%%%%%%%%%%%% layer 2
Q1(i)=-(q1(i)-g.*C(i)).*(hl3(i+1)-hl2(i));
Q2(i)=-(q21(i)-g.*C21(i)).*(h1(i+1)-h1(i));


F21s(i)=(Q1(i)+(B1e(i)));

F42s(i)=(Q2(i)+0.*B21(i)-0.*S1(i));

F22s(i)=(0.*Q1(i)+(B1(i))-0.*Sl1(i));

F41s(i)=(Q2(i)+(B21(i))-S1(i));

FS(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);

FS2(i)=F21(i).*F42(i)-F22(i).*F41(i);
FS41=F;
FS41(i)=FS(i)./(FS2(i));
FS41(isinf(FS41))=0;
FS41(isnan(FS41))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FCE=H;
FCE(i)=F21s(i).*F42s(i)-F22s(i).*F41s(i);

%%%%%%%%%%%%%%%%%
J1(i)=(g./2*(h1(i)./h1(i+1)-h1(i+1)./h1(i)).*(h1(i)-h1(i+1))+0.*g./2.*(h1(i+1)+h1(i))./(h1(i+1).*h1(i)).*(Z1(i)-Z1(i+1)).*(h1(i)-h1(i+1))+...
    -1./(h1(i+1).*h1(i)).*B21e(i).*(h1(i)-h1(i+1)));


ST1(i)=(U(i)./h1(i).*h1(i+1)-sqrt(J1(i)).*h1(i+1));

aj=find(J1<0);

ST1(aj)=0; %#ok<*FNDSB> 


S7(i)=-0.*g./2.*(H(i+2)+H(i-1)).*(Z(i+2)-Z(i-1));

S8(i)=-0.*g./2.*(Hl1(i+2)+Hl1(i-1)).*(Z(i+2)-Z(i-1));

B(i)=0.5.*g.*(H(i-1)+H(i+2)).*rr.*(Hl1(i+2)-Hl1(i-1));

B0(i)=0.5.*g.*(Hl1(i-1)+Hl1(i+2)).*(H(i+2)-H(i-1));

Fp(i)=(0.5.*(U(i-1)+U(i-1))).^2./H(i+2)+0.5.*g.*(H(i+2)).^2;
Fm(i)=(0.5.*(U(i-1)+U(i-1))).^2./H(i-1)+0.5.*g.*(H(i-1)).^2;

Fp2(i)=(0.5.*(Ul1(i-1)+Ul1(i-1))).^2./Hl1(i+2)+0.5.*g.*Hl1(i+2).^2;
Fm2(i)=(0.5.*(Ul1(i-1)+Ul1(i-1))).^2./Hl1(i-1)+0.5.*g.*Hl1(i-1).^2;

k42(i)=-(Fp(i)-Fm(i)-S7(i)+B(i))./(H(i+2)-H(i-1));
k22(i)=(Fp2(i)-Fm2(i)-S8(i)+B0(i))./(H(i+2)-H(i-1));


k21(i)=(Fp2(i)-Fm2(i))./(Hl1(i+2)-Hl1(i-1));

k41(i)=(+B(i))./(Hl1(i+2)-Hl1(i-1));


k(i)=-(k42(i).*k21(i)-k22(i).*k41(i));

aj=find(k<0);
u2(aj)=min(ST1(aj),ST1(aj));
aj=find(k>0);
u2(aj)=max(ST1(aj),ST1(aj));

B21(i)=0.5.*g.*(h1(i+1)+h1(i)).*rr.*(hl01(i+1)-hl01(i));

B21E(i)=0.5.*g.*(h1(i+1)+h1(i)).*rr.*(h1(i)-h1(i+1)+Z1(i)-Z1(i+1));


B2(i)=0.5.*g.*(h3(i+1)+h2(i)).*rr.*(hl03(i+1)-hl02(i));

B1(i)=0.5.*g.*(hl3(i+1)+hl2(i)).*(H(i+1)-H(i));


Fp(i)=u2(i).^2./h1(i+1)+0.5.*g.*(h1(i+1)).^2;

Fm(i)=U(i).^2./h1(i)+0.5.*g.*(h1(i)).^2;

k42(i)=(Fp(i)-Fm(i)+B21e(i));

k22(i)=(Tl11(i)-Sl1(i)+B1(i));

k21(i)=(Tl11(i));

k41(i)=(B21(i));

k(i)=-(k42(i).*k21(i)-k22(i).*k41(i));

k42d(i)=(T2(i)+B2e(i));

k42c(i)=(T1(i)+B21e(i));


DT(i)=(k42d(i)-k42(i))./(k42d(i)-k42c(i));
DT(isnan(DT))=0;
DT(isinf(DT))=0;
DT(DT<0)=0;
DT(DT>1)=1;


F1=F2;
   Tfs(i)=(1-DT(i)).*(T2(i)+B2e(i))+DT(i).*(T1(i)+B21e(i));
F2(i)=(1-DT(i)).*FS3(i)+DT(i).*FS31(i);
F1(i)=(1-DT(i)).*FS4(i)+DT(i).*FS41(i);

F1s=H;
F2s(i)=(1-DT(i)).*FS3(i)+DT(i).*0;
F1s(i)=(1-DT(i)).*FS4(i)+DT(i).*0;


FF(i)=1;
FF(ttY)=0;

G2(i)=Hl1(i+1)-Hl1(i);

G3(i)=0.5.*(1/g).*((0.5.*(U(i+1)+U(i))./H(i+1)).^2-(0.5.*(U(i+1)+U(i))./H(i)).^2)+H(i+1)-H(i)+Z(i+1)-Z(i);

G(i)=Ul1(i+1)-Ul1(i);


FU=H*0+1;

FD(i)=1;
FD(ttC)=0;  
D(i)=1;
D(G2<0 )=0;

M(i)=0;
F(i)=F1(i-1);


G(i)=max(min((H(i+1)-H(i))./(H(i+2)-H(i)),1),0);


M(i)=G(i).*M(i)+(1-G(i)).*F(i);
FV=FU;
  FU(i)=((1-FF(i)).*D(i-1).*M(i)+FF(i).*(1));
 FV(i)=1;

   F2(i)=((1-FF(i)).*D(i-1).*min(F2(i),(0.5.*(M(i)+M(i))))+FF(i).*(FD(i).*F2(i)+(1-FD(i)).*F2s(i)));
   F1(i)=((1-FF(i)).*D(i-1).*min(F1(i),(0.5.*(M(i)+M(i))))+FF(i).*(FD(i).*F1(i)+(1-FD(i)).*F1s(i)));


  F2(i)=(1-FC(i)).*max(F2(i+1),F2(i))+FC(i).*F2(i);

  F1(i)=(1-FC(i)).*max(F1(i+1),F1(i))+FC(i).*F1(i);
%%
%Two layer SWEs fixes
F1(i)=(max(0,min(F1(i),1)));
F2(i)=(max(0,min(F2(i),1)));
FU(i)=(max(0,min(FU(i),1)));

%%
% layer 2

C1(i)=(max(F2(i),F2(i)).*(H(i+1)-H(i)));

%%%%%%%%%%%%%%%%%%%%%%%%% HLL solver %%%%%%%%%%%%%%%%%%%%%%%%

c2(i)=sqrt((H(i)+rr*Hl1(i)).*g); %layer 2

c1(i)=sqrt(Hl1(i).*g); %layer 1

SRe1(i)=max(Ul1(i+1)./Hl1(i+1)+c1(i+1),Ul1(i)./Hl1(i)+c1(i)); %layer 1
SLe1(i)=min(Ul1(i+1)./Hl1(i+1)-c1(i+1),Ul1(i)./Hl1(i)-c1(i));

SRe2(i)=max(U(i+1)./H(i+1)+c2(i+1),U(i)./H(i)+c2(i)); %layer 2
SLe2(i)=min(U(i+1)./H(i+1)-c2(i+1),U(i)./H(i)-c2(i));

SRe(i)=max(0,max(SRe2(i),SRe1(i)));

SLe(i)=min(0,min(SLe2(i),SLe1(i)));

    
        X(i)=(SLe(i).*SRe(i)).*(C1(i));
FHm(i)=(-SRe(i).*(U(i+1)-U(i))+X(i))./(SRe(i)-SLe(i));% HLL layer 2
FHp(i)=(-SLe(i).*(U(i+1)-U(i))+X(i))./(SRe(i)-SLe(i));

X(i)=SLe(i).*SRe(i).*FV(i).*(U(i+1)-U(i));

FUm(i)=(-SRe(i).*(Tfs(i))+X(i))./(SRe(i)-SLe(i)); % HLL layer 2
FUp(i)=(-SLe(i).*(Tfs(i))+X(i))./(SRe(i)-SLe(i));

%%
% layer 2

H(i)=H(i)+dt.*(-1./(dx).*(FHp(i)-FHm(i-1)));
U(i)=U(i)+dt.*(-1./(dx).*(FUp(i)-FUm(i-1)));

%%
C1(i)=max(F1(i),F1(i)).*(Hl1(i+1)-Hl1(i));

X(i)=-(SR(i).*(C1(i)));

FHpl1(i)=(0.*(Ul1(i+1)+Ul1(i))+X(i))./(2);
FHml1(i)=(0.*(Ul1(i+1)+Ul1(i))+X(i))./(2);


X(i)=-SR(i).*(Ul1(i+1)-Ul1(i));

FUpl1(i)=(0+X(i))./(2);
FUml1(i)=(0+X(i))./(2);

Q(i)=Ul1(i+1)-Ul1(i);

L1(i)=(((0.5).*Q(i-1)+(0.5).*Q(i)));


L(i)=0.5.*(Tl11(i)+Tl11(i-1));

S(i)=0.5.*(Sl1(i)+Sl1(i-1));

B(i)=0.5.*(B1(i)+B1(i-1));


   X(i)=(SLe(i).*SRe(i)).*(C1(i));
FHpl1(i)=(-SRe(i).*(Ul1(i+1)-Ul1(i))+X(i))./(SRe(i)-SLe(i));%HLL
FHml1(i)=(-SLe(i).*(Ul1(i+1)-Ul1(i))+X(i))./(SRe(i)-SLe(i));

X(i)=SLe(i).*SRe(i).*FU(i).*(Ul1(i+1)-Ul1(i));

FUpl1(i)=(-SRe(i).*(Tl11(i)+B1e(i))+X(i))./(SRe(i)-SLe(i)); %HLL
FUml1(i)=(-SLe(i).*(Tl11(i)+B1e(i))+X(i))./(SRe(i)-SLe(i));

%%

Hl1(i)=Hl1(i)+dt.*(-1./(dx).*((FHpl1(i)-FHml1(i-1))));
Ul1(i)=Ul1(i)+dt.*(-1./(dx).*((FUpl1(i)-FUml1(i-1))));
%%

% Impose Boundary Condition
H(1) = H(2);
U(1) = U(2);

Hl1(1) = Hl1(2);
Ul1(1) = Ul1(2);

H(n-1:n+2) =H(n-2);             
Hl1(n-1:n+2) = Hl1(n-2);
U(n-1:n+2) = U(n-2);
Ul1(n-1:n+2) = Ul1(n-2);



disp("t=");disp(t)
disp("density ratio=");disp(rr)

end


end
end
