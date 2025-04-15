clc
clear all


h=0.5; %in A
r=48:h:800; %Stern layer 15A plus radius 20A
%for PAH wrapped [Na+ shell size (A)	15	15	15	19	19	19	28]
perm=8.8541878128E-12; %in C/V/m
diele=80;
e=1.602176634E-19; % in C
kB = 1.380649E-23; % in J/K
T =300; % in K 
NA=6.02214076E+23;  % #/mol
%bulk c-Na+ 	0.000448505	0.001445183 0.005431894 0.010415282 0.050282392 0.100116279 0.199784053
%bulk c-Cl-	0.00000000	0.000996678 0.004983389 0.009966777 0.049833887 0.099667774 0.199335548
%net charge at 35A: -75.59	-46.186	-41.94	-38.38287634	-34.795	-29.92015152	-41.02416667
%MD_RDF_cNa0_N=[0.0001308	0.0009669	0.0049715	0.0099689	0.0490882	0.1011327	0.2000605];
%MD_RDF_cCl0_N=[0.0000000	0.0009669	0.0049611	0.0099455	0.0491048	0.1011269	0.2001094];
%PAH_net charge at 35A: -19.5043	-10.6229	-9.1432	-7.0479	-6.5286	-8.3015	-2.8977
%PAH_cNa0[0.0003634357	0.0013948046	0.0053312713	0.0101944438	0.0503253685	0.0988951938	0.1987859474]
%PAH_cCl0[0.0003616542	0.0013967198	0.0053238111	0.0101885543	0.0503214844	0.0989022107	0.1987674982]
Q_N=[-19.5043	-10.6229	-9.1432	 -7.0479	-6.5286	-8.3015	-2.8977];
cNa0_N=[0.0003634357	0.0013948046	0.0053312713	0.0101944438	0.0503253685	0.0988951938	0.1987859474];
cCl0_N=[0.0003616542	0.0013967198	0.0053238111	0.0101885543	0.0503214844	0.0989022107	0.1987674982];
scale=1; % scale Potentail Pr 

for Ni = 1:7    
clear Pr cNa cCl rouNa rouCl
Q=Q_N(1,Ni); 
cNa0 = cNa0_N(1,Ni); % in mol/L
cCl0 = cCl0_N(1,Ni); % in mol/L 

phi0=Q*e/4/pi/perm/diele/(r(1,1)*1E-10); %in V  %BC1 surface potential
phi1=-Q*e*(h*1E-10)/4/pi/perm/diele/(r(1,1)*1E-10)^2+phi0; %in V/m  %BC2 surface E-field 
Pr(1,1)=phi0;
for i = 2:size(r,2)
    Pr(1,i)=Pr(1,i-1)*(phi1/phi0); % initialize the potentials 
end
   %figure; plot(r,Pr(1,:));
%% 1st iterative solving Pr the potential field

rouNa(1,1)=0;
rouCl(1,1)=0;
for i = 2:size(r,2)
    cNa(1,i)=cNa0*exp(-e*Pr(1,i)/scale/(kB*T));
    cCl(1,i)=cCl0*exp(e*Pr(1,i)/scale/(kB*T));
    rouNa(1,i)=cNa0/1000*e*NA*exp(-e*Pr(1,i)/scale/(kB*T));
    rouCl(1,i)=cCl0/1000*e*NA*exp(e*Pr(1,i)/scale/(kB*T));
end

Pr(2,1)=Pr(1,1);Pr(2,2)=Pr(1,2);
for i = 3:size(r,2)
    Pr(2,i)= ((rouNa(1,i)-rouCl(1,i))/perm/diele+Pr(1,i-1)*2/(h*1E-10)^2-Pr(1,i-2)*(1/(h*(1E-10))^2-1/(h*(1E-10)*r(1,i)*(1E-10))))/(1/(h*1E-10)^2+1/(h*(1E-10)*r(1,i)*(1E-10)));
if Pr(2,i)< Pr(2,i-1)||(Pr(2,i)>0)
    Pr(2,i)= Pr(1,i);
end

end

%% more iterations
itn=1500;
for it =2:itn
    
rouNa(it,1)=0;
rouCl(it,1)=0;

for i = 2:size(r,2)
    rouNa(it,i)=cNa0/1000*e*NA*exp(-e*Pr(it,i)/scale/(kB*T));
    rouCl(it,i)=cCl0/1000*e*NA*exp(e*Pr(it,i)/scale/(kB*T));
    cNa(it,i)=cNa0*exp(-e*Pr(it,i)/scale/(kB*T));
    cCl(it,i)=cCl0*exp(e*Pr(it,i)/scale/(kB*T));
end

Pr(it+1,1)=Pr(it,1);Pr(it+1,2)=Pr(it,2);
for i = 3:size(r,2)
    Pr(it+1,i)= ((rouNa(it,i)-rouCl(it,i))/perm/diele+Pr(it,i-1)*2/(h*1E-10)^2-Pr(it,i-2)*(1/(h*(1E-10))^2-1/(h*(1E-10)*r(1,i)*(1E-10))))/(1/(h*1E-10)^2+1/(h*(1E-10)*r(1,i)*(1E-10)));
if (Pr(it+1,i)< Pr(it+1,i-1))||(Pr(it+1,i)>0)
    Pr(it+1,i)= Pr(it,i);
end
end

%if mod(it,1000)==0
 %  figure; plot(r,Pr(it,:));
 %  figure; plot(r,cNa(it,:));
%end
end
    %% a few extra iterations

itnx=100;
for it =itn+1:itn+itnx
    
rouNa(it,1)=0;
rouCl(it,1)=0;

for i = 2:size(r,2)
    rouNa(it,i)=cNa0/1000*e*NA*exp(-e*Pr(it,i)/scale/(kB*T));
    rouCl(it,i)=cCl0/1000*e*NA*exp(e*Pr(it,i)/scale/(kB*T));
    cNa(it,i)=cNa0*exp(-e*Pr(it,i)/scale/(kB*T));
    cCl(it,i)=cCl0*exp(e*Pr(it,i)/scale/(kB*T));
end

Pr(it+1,1)=Pr(it,1);Pr(it+1,2)=Pr(it,2);
for i = 3:size(r,2)
    Pr(it+1,i)= ((rouNa(it,i)-rouCl(it,i))/perm/diele+Pr(it+1,i-1)*2/(h*1E-10)^2-Pr(it+1,i-2)*(1/(h*(1E-10))^2-1/(h*(1E-10)*r(1,i)*(1E-10))))/(1/(h*1E-10)^2+1/(h*(1E-10)*r(1,i)*(1E-10)));
%if Pr(it+1,i)>0
 %   Pr(it+1,i)= 0;
%end
end

end
 resultPr(:,Ni)=Pr(:,931);
 result(:,1)=transpose(r);
 result(:,Ni+1)=transpose(Pr(itn+itnx+1,:));
 resultNa(:,1)=transpose(r);
 resultNa(:,Ni+1)=transpose(cNa(it,:));
 resultCl(:,1)=transpose(r);
 resultCl(:,Ni+1)=transpose(cCl(it,:));
end
   figure; plot(resultPr);
   figure; plot(result(:,1),result(:,2:8));
   figure; plot(resultNa(:,1),resultNa(:,2:8));
   figure; plot(resultCl(:,1),resultCl(:,2:8));