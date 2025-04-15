clear all;
clc

%90 citrate3- on the Au surface
filename = 'Na.txt';%change data file name;
B = importdata(filename);

com=[48.5474 48.5474 48.5474];%com of AuNP

nt=50;
ni=size(B,1)/nt;

for i=1:nt
    t(1,i)=i*20000/1000000; %ns
    for j = 1:ni
        xt=B(j+ni*(i-1),1);
        yt=B(j+ni*(i-1),2);
        zt=B(j+ni*(i-1),3);
        molt(j,1)=xt;molt(j,2)=yt;molt(j,3)=zt;
        posm(j,i)=sqrt((xt-com(1,1))^2+(yt-com(1,2))^2+(zt-com(1,3))^2);
    end
end

%RDF 
clear r kc nc
nn=1600; %rdf bin numbers
Rmax=800; %rdf edge
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
dr=Rmax/nn; r1=10; r2=50;r3=100;r4=200;r5=300;r6=400;
r=dr:dr:Rmax;
for i = 1:nn
    coef(i,1)=4*pi*r(1,i)^2*dr;
    rdf(i,:)=nc(i,:)/coef(i,1);
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
set(plot1(1),'DisplayName','50 Å');
set(plot1(2),'DisplayName','100 Å');
set(plot1(3),'DisplayName','200 Å');
set(plot1(4),'DisplayName','300 Å');
set(plot1(5),'DisplayName','400 Å');

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
    aveRDF_Na(i,1)=mean(rdf(i,1:nt));
end
subplot(3,1,2)
plot(r,aveRDF_Na,'r'); %last 10 frames rdf,r,rdf(:,nt-6:nt),
xlabel('r(Å)','FontSize',11);
ylabel('g(r)','FontSize',11);

for i=1:nn
    aveNum_Na(i,1)=mean(nc(i,1:nt));
    sumNum_Na(i,1)=sum(aveNum_Na(1:i,1));
end
subplot(3,1,3)
plot(r,sumNum_Na,'r'); %r,nc(:,nt-6:nt),
xlabel('r(Å)','FontSize',11); %last 10 frames number count
ylabel('Na+ number','FontSize',11);
%%
result(:,1)=transpose(r);result(:,2)=sumNum_Na;
%%
clear ni rdf;

filename2 = 'Cl.txt';%change data file name;
B2 = importdata(filename2);
ni=size(B2,1)/nt;
for i=1:nt
    t(1,i)=i*20000/1000000; %ns
    for j = 1:ni
        xt=B2(j+ni*(i-1),1);
        yt=B2(j+ni*(i-1),2);
        zt=B2(j+ni*(i-1),3);
        molt(j,1)=xt;molt(j,2)=yt;molt(j,3)=zt;
        posm2(j,i)=sqrt((xt-com(1,1))^2+(yt-com(1,2))^2+(zt-com(1,3))^2);
    end
end

% RDF 
clear r kc nc
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
dr=Rmax/nn; r1=10; r2=50;r3=100;r4=200;r5=300;r6=400;
r=dr:dr:Rmax;
for i = 1:nn
    coef(i,1)=4*pi*r(1,i)^2*dr;
    rdf(i,:)=nc2(i,:)/coef(i,1);
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
set(plot1(1),'DisplayName','50 Å');
set(plot1(2),'DisplayName','100 Å');
set(plot1(3),'DisplayName','200 Å');
set(plot1(4),'DisplayName','300 Å');
set(plot1(5),'DisplayName','400 Å');

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
    nc0(1,i)=kc;
end
%
figure;
subplot(3,1,1)
hist(posm2(:,nt),50); %last frame histogram
xlabel('r(Å)','FontSize',11);
ylabel('Cl- number','FontSize',11);

hold on;

for i=1:nn
    aveRDF_Cl(i,1)=mean(rdf(i,1:nt));
end
subplot(3,1,2)
plot(r,aveRDF_Cl,'r'); %last 10 frames rdf,r,rdf(:,nt-6:nt),
xlabel('r(Å)','FontSize',11);
ylabel('g(r)','FontSize',11);

for i=1:nn
    aveNum_Cl(i,1)=mean(nc2(i,1:nt));
    sumNum_Cl(i,1)=sum(aveNum_Cl(1:i,1));
end
subplot(3,1,3)
plot(r,sumNum_Cl,'r'); %r,nc(:,nt-6:nt),
xlabel('r(Å)','FontSize',11); %last 10 frames number count
ylabel('Cl- number','FontSize',11);
%%
figure;
plot(r, sumNum_Na-sumNum_Cl, 'r');
result(:,3)=sumNum_Cl; result(:,4)=sumNum_Na-sumNum_Cl;