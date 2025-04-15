clear all;
clc

%90 citrate3- on the Au surface
filename = 'Na.txt';%change data file name;
B = importdata(filename);

filename1 = 'log.8pah0.001M_nvt22';%change data file name;
A = importdata(filename1, ' ', 11989); %406 %1459 %2512 %3565 %4618 %5671 %6724 %7777 %8830 %9883  %10936 %11989 %13042 %14095 %15148 %16201 %17254 
C=A.data;

nt=101;
ni=size(B,1)/nt;
dn = 10000; %output steps
ts = 1000; %timestep 1000
com(1,1:3)=C(1,7:9);%com of AuNP
for i=2:nt
    ii = (i-1)*dn/ts;
com(i,1:3)=C(ii,7:9);%com of AuNP
end

for i=1:nt
    t(1,i)=i*dn/1000000; %ns
    for j = 1:ni
        xt=B(j+ni*(i-1),1);
        yt=B(j+ni*(i-1),2);
        zt=B(j+ni*(i-1),3);
        molt(j,1)=xt;molt(j,2)=yt;molt(j,3)=zt;
        posm(j,i)=sqrt((xt-com(i,1))^2+(yt-com(i,2))^2+(zt-com(i,3))^2);
    end
end

%RDF 
clear r kc nc
nn=1600; %rdf bin numbers
Rmax=800; %rdf edge
NA=6.0221408E+23;
for i = 1:nt
    kc(1:nn,1)=0;
    for j = 1:ni
        temp = posm(j,i);
        for ii = 1:nn
        if temp < ii*(Rmax/nn)
            kc(ii,1)=kc(ii,1)+1;
            break;
        end
        end
    end
    nc(:,i)=kc;
end
%
dr=Rmax/nn; r1=10; r2=100;r3=150;r4=200;r5=250;r6=300;
r=dr:dr:Rmax;
for i = 1:nn
    coef(i,1)=4*pi*r(1,i)^2*dr;
    rdf_na(i,:)=nc(i,:)/coef(i,1)*1e27/NA; %unit convert to M
end
n1=round(r1/dr);
n2=round(r2/dr);n3=round(r3/dr);n4=round(r4/dr);n5=round(r5/dr);n6=round(r6/dr);
for i = 1:size(posm,2)
    Na(1,i)=sum(nc(n1:n2,i));
    Na(2,i)=sum(nc(n1:n3,i));
    Na(3,i)=sum(nc(n1:n4,i));
    Na(4,i)=sum(nc(n1:n5,i));
    Na(5,i)=sum(nc(n1:n6,i));
end
%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(t,Na,'Parent',axes1);
set(plot1(1),'DisplayName','R = 100 Å');
set(plot1(2),'DisplayName','R = 150 Å');
set(plot1(3),'DisplayName','R = 200 Å');
set(plot1(4),'DisplayName','R = 250 Å');
set(plot1(5),'DisplayName','R = 300 Å');

box(axes1,'on');
legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5 0.5 0 0]);
xlabel('Simulation time (ns)','FontSize',11); %last 10 frames number count
ylabel('Total Na+ number in the shell','FontSize',11);

%
radm=abs(posm-20);
for i = 1:size(radm,2)
    kc=0;
    for j = 1:size(radm,1)
        if radm(j,i)<75
            kc=kc+1;
        end
    end
    nc0(1,i)=kc;
end
%
figure;
%plot(r,nc(:,t1));
%subplot(3,1,2)
%plot(r,rdf(:,5));
%plot(r,nc(:,t2));
%subplot(3,1,3)
%plot(t,Na);
subplot(3,1,1)
hist(posm(:,nt),50); %last frame histogram
xlabel('r(Å)','FontSize',11);
ylabel('Na+ number','FontSize',11);

hold on;

for i=1:nn
    aveRDF_Na(i,1)=mean(rdf_na(i,2:nt));
end
subplot(3,1,2)
plot(r,aveRDF_Na,'r'); %last 10 frames rdf,r,rdf(:,nt-6:nt),
xlabel('r(Å)','FontSize',11);
ylabel('g(r)','FontSize',11);

for i=1:nn
    aveNum_Na(i,1)=mean(nc(i,2:nt));
    sumNum_Na(i,1)=sum(aveNum_Na(1:i,1));
end
subplot(3,1,3)
plot(r,sumNum_Na,'r'); %r,nc(:,nt-6:nt),
xlabel('r(Å)','FontSize',11); %last 10 frames number count
ylabel('Na+ number','FontSize',11);
%%
clear ni;

filename2 = 'Cl.txt';%change data file name;
B2 = importdata(filename2);
ni=size(B2,1)/nt;
for i=1:nt
    t(1,i)=i*dn/1000000; %ns
    for j = 1:ni
        xt=B2(j+ni*(i-1),1);
        yt=B2(j+ni*(i-1),2);
        zt=B2(j+ni*(i-1),3);
        molt(j,1)=xt;molt(j,2)=yt;molt(j,3)=zt;
        posm2(j,i)=sqrt((xt-com(1,1))^2+(yt-com(1,2))^2+(zt-com(1,3))^2);
    end
end

% RDF 
clear r kc;
nn=1600; %rdf bin numbers
Rmax=800; %rdf edge
for i = 1:nt
    kc(1:nn,1)=0;
    for j = 1:ni
        temp = posm2(j,i);
        for ii = 1:nn
        if temp < ii*(Rmax/nn)
            kc(ii,1)=kc(ii,1)+1;
            break;
        end
        end
    end
    nc2(:,i)=kc;
end
%% 
dr=Rmax/nn; r1=10; r2=100;r3=150;r4=200;r5=250;r6=300;
r=dr:dr:Rmax;
for i = 1:nn
    coef(i,1)=4*pi*r(1,i)^2*dr; %unit convert to M unit
    rdf_cl(i,:)=nc2(i,:)/coef(i,1)*1e27/NA;
end
n1=round(r1/dr);
n2=round(r2/dr);n3=round(r3/dr);n4=round(r4/dr);n5=round(r5/dr);n6=round(r6/dr);
for i = 1:size(posm,2)
    Cl(1,i)=sum(nc2(n1:n2,i));
    Cl(2,i)=sum(nc2(n1:n3,i));
    Cl(3,i)=sum(nc2(n1:n4,i));
    Cl(4,i)=sum(nc2(n1:n5,i));
    Cl(5,i)=sum(nc2(n1:n6,i));
end
%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(t,Cl,'Parent',axes1);
set(plot1(1),'DisplayName','100 Å');
set(plot1(2),'DisplayName','150 Å');
set(plot1(3),'DisplayName','200 Å');
set(plot1(4),'DisplayName','250 Å');
set(plot1(5),'DisplayName','300 Å');

box(axes1,'on');
legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5 0.5 0 0]);
xlabel('Simulation time (ns)','FontSize',11); %last 10 frames number count
ylabel('Total Cl- number in the shell','FontSize',11);
%
radm=abs(posm2-20);
for i = 1:size(radm,2)
    kc=0;
    for j = 1:size(radm,1)
        if radm(j,i)<75
            kc=kc+1;
        end
    end
    nc0(2,i)=kc;
end
%
figure;
subplot(3,1,1)
hist(posm2(:,nt),50); %last frame histogram
xlabel('r(Å)','FontSize',11);
ylabel('Cl- number','FontSize',11);

hold on;

for i=1:nn
    aveRDF_Cl(i,1)=mean(rdf_cl(i,2:nt));
end
subplot(3,1,2)
plot(r,aveRDF_Cl,'r'); %last 10 frames rdf,r,rdf(:,nt-6:nt),
xlabel('r(Å)','FontSize',11);
ylabel('g(r)','FontSize',11);

for i=1:nn
    aveNum_Cl(i,1)=mean(nc2(i,2:nt));
    sumNum_Cl(i,1)=sum(aveNum_Cl(1:i,1));
end
subplot(3,1,3)
plot(r,sumNum_Cl,'r'); %r,nc(:,nt-6:nt),
xlabel('r(Å)','FontSize',11); %last 10 frames number count
ylabel('Cl- number','FontSize',11);


%%
clear ni;

filename3 = 'N.txt';%change data file name;
B3 = importdata(filename3);
ni=size(B3,1)/nt;
for i=1:nt
    t(1,i)=i*dn/1000000; %ns
    for j = 1:ni
        xt=B3(j+ni*(i-1),1);
        yt=B3(j+ni*(i-1),2);
        zt=B3(j+ni*(i-1),3);
        molt(j,1)=xt;molt(j,2)=yt;molt(j,3)=zt;
        posm3(j,i)=sqrt((xt-com(1,1))^2+(yt-com(1,2))^2+(zt-com(1,3))^2);
    end
end

% RDF 
clear r kc;
nn=1600; %rdf bin numbers
Rmax=800; %rdf edge
for i = 1:nt
    kc(1:nn,1)=0;
    for j = 1:ni
        temp = posm3(j,i);
        for ii = 1:nn
        if temp < ii*(Rmax/nn)
            kc(ii,1)=kc(ii,1)+1;
            break;
        end
        end
    end
    nc3(:,i)=kc;
end
%
dr=Rmax/nn; r1=10; r2=100;r3=150;r4=200;r5=250;r6=300;
r=dr:dr:Rmax;
for i = 1:nn
    coef(i,1)=4*pi*r(1,i)^2*dr; %unit convert to M unit
    rdf_nh3(i,:)=nc3(i,:)/coef(i,1)*1e27/NA;
end
n1=round(r1/dr);
n2=round(r2/dr);n3=round(r3/dr);n4=round(r4/dr);n5=round(r5/dr);n6=round(r6/dr);
for i = 1:size(posm,2)
    N(1,i)=sum(nc3(n1:n2,i));
    N(2,i)=sum(nc3(n1:n3,i));
    N(3,i)=sum(nc3(n1:n4,i));
    N(4,i)=sum(nc3(n1:n5,i));
    N(5,i)=sum(nc3(n1:n6,i));
end
%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(t,N,'Parent',axes1);
set(plot1(1),'DisplayName','100 Å');
set(plot1(2),'DisplayName','150 Å');
set(plot1(3),'DisplayName','200 Å');
set(plot1(4),'DisplayName','250 Å');
set(plot1(5),'DisplayName','300 Å');

box(axes1,'on');
legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5 0.5 0 0]);
xlabel('Simulation time (ns)','FontSize',11); %last 10 frames number count
ylabel('Total Cl- number in the shell','FontSize',11);
%
radm=abs(posm3-20);
for i = 1:size(radm,2)
    kc=0;
    for j = 1:size(radm,1)
        if radm(j,i)<75
            kc=kc+1;
        end
    end
    nc0(3,i)=kc;
end
%
figure;
subplot(3,1,1)
hist(posm3(:,nt),50); %last frame histogram
xlabel('r(Å)','FontSize',11);
ylabel('NH3+ number','FontSize',11);

hold on;

for i=1:nn
    aveRDF_N(i,1)=mean(rdf_nh3(i,2:nt));
end
subplot(3,1,2)
plot(r,aveRDF_N,'r'); %last 10 frames rdf,r,rdf(:,nt-6:nt),
xlabel('r(Å)','FontSize',11);
ylabel('g(r)','FontSize',11);

for i=1:nn
    aveNum_N(i,1)=mean(nc3(i,2:nt));
    sumNum_N(i,1)=sum(aveNum_N(1:i,1));
end
subplot(3,1,3)
plot(r,sumNum_N,'r'); %r,nc(:,nt-6:nt),
xlabel('r(Å)','FontSize',11); %last 10 frames number count
ylabel('NH3+ number','FontSize',11);

%% E field and potential
Qe = sumNum_Na-sumNum_Cl+sumNum_N-270;
perm0=8.8541878128E-12; %vacuum permeability in(F/m) or (C/V/m) or (C^2/N/m^2)
perm1=80; %relative permeability of water 
e0 =1.602176634E-19; %electron charge in C
for i = 1:nn
    Ef(i,1)=Qe(i,1)*e0/(4*pi*perm0*perm1*(r(1,i)/1e10)^2);
end

Ef(nn-i+1,2)=0;
for i = 1:nn-1
    Ef(nn-i,2)=dr/1e10*0.5*(Ef(nn-i,1)+Ef(nn-i+1,1))+Ef(nn-i+1,2);
end

%%
figure;
subplot(2,2,1)
plot(r, Qe, 'r');
xlim([20,500]);
xlabel('r(Å)'); ylabel('Total net charge Q(r)');
subplot(2,2,2)
plot(r, aveRDF_Na);
xlim([20,500]);ylim([0,0.005]);
xlabel('r(Å)'); ylabel('Na+ density (M)');
subplot(2,2,3)
plot(r, aveRDF_Cl);
xlim([20,500]);ylim([0,0.005]);
xlabel('r(Å)'); ylabel('Cl- density (M)');
subplot(2,2,4)
plot(r, aveRDF_N);
xlim([20,500]);%ylim([0,0.005]);
xlabel('r(Å)'); ylabel('NH3+ density (M)');

figure;
subplot(2,1,1)
plot(r, Ef(:,1), 'r');
xlim([20,500]);
xlabel('r(Å)'); ylabel('E(r) electric field V/m');
subplot(2,1,2)
plot(r, Ef(:,2));
xlim([20,500]);
xlabel('r(Å)'); ylabel('phi(r) electric potential V');

result(:,1)=transpose(r);
result(:,2)=aveRDF_Na;
result(:,3)=aveRDF_Cl; 
result(:,4)=sumNum_Na;
result(:,5)=sumNum_Cl;
result(:,6)=Qe;
result(:,7)=Ef(:,1);
result(:,8)=Ef(:,2);

%transpose(Na)
transpose(Cl)