
%%%%%%%%%%%
% Forward model for generating idealized dataset of 228,230,234Th and
% Particles. Model equations are discretized using a backwards differences
% scheme. This routine also generate errors for the idealized data,
% assuming the relative error in the data is 10% for 228,230,234Th, 228Ra, and particle
% data, and 5% for 234,238U. Finally this routine also generates noise
% contaminated idealized dataset, which are used as a starting point for
% the inversion performed by FMINCON in "mainThP.m"  

% Last modified by Paul Lerner
% 02/13/2017
%%%%%%%%%%%




clear all
close all;
clc;

% depth from 100m, every 100 m, to 4100 m.
z=[100:100:4100];

% decay constants
dTh230=log(2)/75690;
dTh228=log(2)/1.91;
dTh234=log(2)/(24.1/365.242);

% rate constant values for adsorption (k1), desorption (kdes), particle
% degradation (Brem), and sinking speed (w).
kdes=3;
k1=0.7;
Brem=1;
w=600;

% boundary values for Th and particles
Th234p(1)=600;
Th230p(1)=0.0005;
Th228p(1)=0.8;
P(1)=15;

% model assumes constatn 234,238U
U238=2400;
U234=2750;

% calculate Ra228 as a linear combination of two exponentials
A=30;
B=5;
Ra228=(A.*exp(-(z-z(1))/500))+B.*exp(-(z(end)-z)/500);

% depth difference is 100 m
dz=100;

% initialize x, for which we attempt to find a solution. x will contain the
% Th isotope activities and particle ocncentrations
x=zeros(length(z),7);
x(1,2)=Th234p(1);
x(1,4)=Th230p(1);
x(1,6)=Th228p(1);
x(1,7)=P(1);


% Ax=b, where A is the design matrix that describes the model, b are the
% forcings, and x are the unknown Th isotope activites and particle
% concentrations for which we seek a solution. We do this iteratively, updating the forcing vector which contains
% the particulate values at the previous iteration (first iteration uses the boundary value) at every step.
A=[-(k1+dTh234) (kdes+Brem) 0 0 0 0 0;k1 -(kdes+Brem+dTh234+(w/dz)) 0 0 0 0 0;0 0 -(k1+dTh230) (kdes+Brem) 0 0 0;0 0 k1 -(kdes+Brem+dTh230+(w/dz)) 0 0 0;0 0 0 0 -(k1+dTh228) (kdes+Brem) 0;0 0 0 0 k1 -(kdes+Brem+dTh228+(w/dz)) 0;0 0 0 0 0 0 -(Brem+(w/dz))];
b=[-dTh234*U238;-Th234p(1)*w/dz;-dTh230*U234;-Th230p(1)*w/dz;-dTh228*Ra228(1);-Th228p(1)*w/dz;-P(1)*w/dz];

for i=2:length(z);
b(5)=-dTh228*Ra228(i);    
x(i,1:7)=A\b;

b(2)=-x(i,2)*w/dz;
b(4)=-x(i,4)*w/dz;
b(6)=-x(i,6)*w/dz;
b(7)=-x(i,7)*w/dz;
end

% extract each Th isotope and particle concentration from x, and define
% vector of 234,238U.
U238model=U238.*ones(length(z)-1,1);
U234model=U234.*ones(length(z)-1,1);
Ra228model=Ra228(2:end);
Th234dmodel=x(2:end,1);
Th234pmodel=x(:,2);
Th230dmodel=x(2:end,3);
Th230pmodel=x(:,4);
Th228dmodel=x(2:end,5);
Th228pmodel=x(:,6);
Pmodel=x(:,7);


% relative errors in the idealized data are 10%. Also, this relative error
% is used to generate noise contaminated data. 
for i=1:length(U238model);
U238modelerr(i)=U238model(i)*.05;
U238modelrnd(i)=U238model(i)+(U238model(i)/20.*randn(1));
end
for i=1:length(U234model);
U234modelerr(i)=U234model(i)*.05;
U234modelrnd(i)=U234model(i)+(U234model(i)/20.*randn(1));
end
for i=1:length(Ra228model);
Ra228modelerr(i)=Ra228model(i)*.1;
Ra228modelrnd(i)=Ra228model(i)+(Ra228model(i)/10.*randn(1));
end
for i=1:length(Th234dmodel);
Th234dmodelerr(i)=Th234dmodel(i)*.1;
Th234dmodelrnd(i)=Th234dmodel(i)+(Th234dmodel(i)/10.*randn(1));
end
for i=1:length(Th234pmodel);
Th234pmodelerr(i)=Th234pmodel(i)*.1;
Th234pmodelrnd(i)=Th234pmodel(i)+(Th234pmodel(i)/10.*randn(1));
end

for i=1:length(Th230dmodel);
Th230dmodelerr(i)=Th230dmodel(i)*.1;
Th230dmodelrnd(i)=Th230dmodel(i)+(Th230dmodel(i)/10.*randn(1));
end
for i=1:length(Th230pmodel);
Th230pmodelerr(i)=Th230pmodel(i)*.1;
Th230pmodelrnd(i)=Th230pmodel(i)+(Th230pmodel(i)/10.*randn(1));
end
for i=1:length(Th228dmodel);
Th228dmodelerr(i)=Th228dmodel(i)*.1;
Th228dmodelrnd(i)=Th228dmodel(i)+(Th228dmodel(i)/10.*randn(1));
end
for i=1:length(Th228pmodel);
Th228pmodelerr(i)=Th228pmodel(i)*.1;
Th228pmodelrnd(i)=Th228pmodel(i)+(Th228pmodel(i)/10.*randn(1));
end
for i=1:length(Pmodel);
Pmodelerr(i)=Pmodel(i)*.1;
Pmodelrnd(i)=Pmodel(i)+(Pmodel(i)/10.*randn(1));
end

save('Thmodel.mat','U238model','U234model','Ra228model','Th234dmodel','Th234pmodel','Th230dmodel','Th230pmodel','Th228dmodel','Th228pmodel','Pmodel'...
    ,'U238modelerr','U234modelerr','Ra228modelerr','Th234dmodelerr','Th234pmodelerr','Th230dmodelerr','Th230pmodelerr','Th228dmodelerr','Th228pmodelerr','Pmodelerr','z'...
    ,'U238modelrnd','U234modelrnd','Ra228modelrnd','Th234dmodelrnd','Th234pmodelrnd','Th230dmodelrnd','Th230pmodelrnd','Th228dmodelrnd','Th228pmodelrnd','Pmodelrnd');