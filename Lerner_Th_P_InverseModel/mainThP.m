%%%%%%%%%%%
% Algorithm for estimating the rate constants for thorium and particle
% cycling. Fits a single particle class model for thorium and particle
% cycling to 228Th, 230Th, and/or 234Th, their resepctive radioactive parents,
% and particle data using a constrained
% optimization routine (either the Algorithm of Total Inversion or FMINCON).
% The output will be the rate constants for thorium adsorption
% (k1), thorium desorption (k-1), particle degredation (B-1), and particle
% sinking speed (w), as well as poseterior estimates of the thorium isotope
% activities, their parent isotope acitivites, and particle concentration.

% For combinations:
% 1: 234Th
% 2: 230Th
% 3: 228Th
% 4: 234,230Th
% 5: 234,228Th
% 6: 228,230Th
% 7: 228,230,234Th

% Last modified by Paul Lerner
% 02/13/2017
%%%%%%%%%%%
clear all
close all
% First, load the data
fprintf('Please load the data\n')
uiopen

global Th234davg Th234pavg U238avg Th230davg Th230pavg U234avg Th228davg Th228pavg Ra228avg Pavg...
    Th234dnorm Th234pnorm U238norm Th230dnorm Th230pnorm U234norm Th228dnorm Th228pnorm Ra228norm Pnorm...
    Th234derrnorm Th234perrnorm U238errnorm Th230derrnorm Th230perrnorm U234errnorm Th228derrnorm Th228perrnorm Ra228errnorm Perrnorm...
    cmbTh Pdepth Ddepth Thdq Thdq2 dTh230 dTh228 dTh234 x0 C0 sumTh234d sumTh234p sumTh230d sumTh230p sumTh228d sumTh228p sumP  Pdo;

% decay constants here
dTh230=log(2)/75690;
dTh228=log(2)/1.91;
dTh234=log(2)/(24.1/365.242);

% Next, choose the combination of Th isotope you are using

x = inputdlg('Enter number corresponding to combination of Th isotopes used. See documentation for which numbers correspond to which isotope combinations',...
             'Combination', [1 50]);
         
cmbTh = str2num(x{:});

x = inputdlg('Enter 1 if the boundary value is the first depth of your dataset for dissolved and particulate data and 0 if your particulate dataset includes a boundary value that is not in your dissolved dataset',...
    'Input', [1 50]);
         
setboundary = str2num(x{:});

x = inputdlg('Enter 0 if you have error covariance matrices and 1 if you have no error covariance (errors stored in array)',...
    'Input', [1 50]);
         
nocovmat = str2num(x{:});

x = inputdlg('Enter 0 if you want to set your initital estimate of x to x0, or 1 if you want to select a different initial estimate',...
    'Input', [1 50]);
         
setini = str2num(x{:})

% Include particles?

x=inputdlg('Enter 0 if you do not have particle concentration data, and 1 if you do',...
    'Input', [1 50]);

Pdo=str2num(x{:});

% Next, select your variables; both for the data vector x0 and the initial
% estimate of the state vector, xk0.

if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234d=input('Please input dissolved Th234 for x0 \n');
Th234p=input('Please input particulate Th234 for x0 \n');
U238=input('Please input U238 for x0  \n');

Th234derrmat=input('Please input dissolved Th234 errors or error covariance matrix \n');
Th234perrmat=input('Please input particulate Th234 errors or error covariance matrix \n');
U238errmat=input('Please input U238 errors or error covariance matrix \n');
end

if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230d=input('Please input dissolved Th230 for x0 \n');
Th230p=input('Please input particulate Th230 for x0 \n');
U234=input('Please input U234 for x0 \n');

Th230derrmat=input('Please input dissolved Th230 errors or error covariance matrix \n');
Th230perrmat=input('Please input particulate Th230 errors or error covariance matrix \n');
U234errmat=input('Please input U234 errors or error covariance matrix \n');
end

if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228d=input('Please input dissolved Th228 for x0 \n');
Th228p=input('Please input particulate Th228 for x0 \n');
Ra228=input('Please input Ra228 for x0 \n');

Th228derrmat=input('Please input dissolved Th228 errors or error covariance matrix \n');
Th228perrmat=input('Please input particulate Th228 errors or error covariance matrix \n');
Ra228errmat=input('Please input Ra228 errors or error covariance \n');
end

if Pdo==1;
P=input('Please inpute particle concentration for x0 \n');
Perrmat=input('Please inpute particle concentration errors or error covariance matrix \n');
end
% grid of dissolved data. Your particulate depth vector should include one more element than your dissolved depth vector, which will be the depth at the upper boundary layer;
Thdq=input('Please inpute vector of particulate thorium/particle concentration depths \n');
Thdq2=input('Please inpute vector of dissolved thorium depths \n');
Ddepth=length(Thdq2);
Pdepth=length(Thdq);
%%

%if setting a different initial estimate of state vector than x0, input
%those values here

if setini==1;
if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234dini=input('Please input dissolved Th234 for xk0 \n');
Th234pini=input('Please input particulate Th234 for xk0 \n');
U238ini=input('Please input U238 for xk0  \n');
end

if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230dini=input('Please input dissolved Th230 for xk0 \n');
Th230pini=input('Please input particulate Th230 for xk0 \n');
U234ini=input('Please input U234 for xk0 \n');
end

if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228dini=input('Please input dissolved Th228 for xk0 \n');
Th228pini=input('Please input particulate Th228 for xk0 \n');
Ra228ini=input('Please input Ra228 for xk0 \n');
end
if Pdo==1;

Pini=input('Please inpute particle concentration for xk0 \n');
end
end
%% 

% if there are no error covariances, diagonalize the error arrays
if nocovmat==1;
if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7    
U238errmat=diag(U238errmat.^2);
Th234derrmat=diag(Th234derrmat.^2);
Th234perrmat=diag(Th234perrmat.^2);
end
if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
U234errmat=diag(U234errmat.^2);    
Th230derrmat=diag(Th230derrmat.^2);
Th230perrmat=diag(Th230perrmat.^2);
end
if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Ra228errmat=diag(Ra228errmat.^2);    
Th228derrmat=diag(Th228derrmat.^2);
Th228perrmat=diag(Th228perrmat.^2);
end
if Pdo==1;
Perrmat=diag(Perrmat.^2);
end
end
%%

% since the first value is the boundary value, dissolved data should start
% at second value in array.
if setboundary==1;


if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234d=Th234d(2:end);
U238=U238(2:end);
Th234derrmat=Th234derrmat(2:end,2:end);
U238errmat=U238errmat(2:end,2:end);
end

if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230d=Th230d(2:end);
U234=U234(2:end);
Th230derrmat=Th230derrmat(2:end,2:end);
U234errmat=U234errmat(2:end,2:end);
end


if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228d=Th228d(2:end);
Ra228=Ra228(2:end);
Th228derrmat=Th228derrmat(2:end,2:end);
Ra228errmat=Ra228errmat(2:end,2:end);
end

if setini==1;
if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234dini=Th234dini(2:end);
U238ini=U238ini(2:end);
end

if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230dini=Th230dini(2:end);
U234ini=U234ini(2:end);
end

if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228dini=Th228dini(2:end);
Ra228ini=Ra228ini(2:end);
end
end
end

%%

% make sure all data are column vector

if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234d=Th234d(:);
Th234p=Th234p(:);
U238=U238(:);
end
if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230d=Th230d(:);
Th230p=Th230p(:);
U234=U234(:);
end
if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228d=Th228d(:);
Th228p=Th228p(:);
Ra228=Ra228(:);
end
if Pdo==1;
P=P(:);
end
if setini==1;
if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234dini=Th234dini(:);
Th234pini=Th234pini(:);
U238ini=U238ini(:);
end
if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230dini=Th230dini(:);
Th230pini=Th230pini(:);
U234ini=U234ini(:);
end
if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228dini=Th228dini(:);
Th228pini=Th228pini(:);
Ra228ini=Ra228ini(:);
end
if Pdo==1;
Pini=Pini(:);
end
end
    

%% normalize activities and concentrations to average values 


    
if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
Th234davg=mean(Th234d);Th234dnorm=Th234d./Th234davg;
Th234pavg=mean(Th234p);Th234pnorm=Th234p./Th234pavg;
U238avg=mean(U238);U238norm=U238./U238avg;

Th234derrnorm=Th234derrmat./(Th234davg.^2);
Th234perrnorm=Th234perrmat./(Th234pavg.^2);
U238errnorm=U238errmat./(U238avg.^2);
end

if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
Th230davg=mean(Th230d);Th230dnorm=Th230d./Th230davg;
Th230pavg=mean(Th230p);Th230pnorm=Th230p./Th230pavg;
U234avg=mean(U234);U234norm=U234./U234avg;

Th230derrnorm=Th230derrmat./(Th230davg.^2);
Th230perrnorm=Th230perrmat./(Th230pavg.^2);
U234errnorm=U234errmat./(U234avg.^2);
end

if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
Th228davg=mean(Th228d);Th228dnorm=Th228d./Th228davg;
Th228pavg=mean(Th228p);Th228pnorm=Th228p./Th228pavg;
Ra228avg=mean(Ra228);Ra228norm=Ra228./Ra228avg;

Th228derrnorm=Th228derrmat./(Th228davg.^2);
Th228perrnorm=Th228perrmat./(Th228pavg.^2);
Ra228errnorm=Ra228errmat./(Ra228avg.^2);
end

if Pdo==1;
Pavg=mean(P);Pnorm=P./Pavg;
Perrnorm=Perrmat./(Pavg.^2);
end

% construct prior estimates of state vector
priorstatevec;


    

%% model equation normalization constants.

[sumTh234d,sumTh234p,sumTh230d,sumTh230p,sumTh228d,sumTh228p,sumP]=modelnorm(x0);

%% Here is wher we perform the inversion

% lower bound for x
lb=10^-6.*ones(length(x0),1);

options=optimoptions('fmincon','Display','iter','MaxIterations',10000,'MaxFunctionEvaluations',4000000,'FunctionTolerance',10^-3,'ConstraintTolerance',10^-4,'SpecifyObjectiveGradient',true);
% choose an "starting point" for the state vector x, xk0. This os x0 if setini is 0. 
if Pdo~=1;
    Pini=[];
end
if setini==1;
 if cmbTh==1
     xk0=[U238ini;Th234dini;Th234pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;

 end
 if cmbTh==2
     xk0=[U234ini;Th230dini;Th230pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;

end
 if cmbTh==3
     xk0=[Ra228ini;Th228dini;Th228pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;


end
 if cmbTh==4
     xk0=[U238ini;U234ini;Th234dini;Th234pini;Th230dini;Th230pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;

 end
 if cmbTh==5
     xk0=[U238ini;Ra228ini;Th234dini;Th234pini;Th228dini;Th228pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;

 end
 if cmbTh==6
     xk0=[U234ini;Ra228ini;Th230dini;Th230pini;Th228dini;Th228pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;

 end
 if cmbTh==7
     xk0=[U238ini;U234ini;Ra228ini;Th234dini;Th234pini;Th230dini;Th230pini;Th228dini;Th228pini;Pini;k1mean;kdesmean;Bremmean;wmean]./x0mean;


 end
else
xk0=x0;
end

[x,fval,exitflag,output]=fmincon('costfmincon',xk0,[],[],[],[],lb,[],'Thpmodel',options);

%% Calculate the error covariance matrix for x

F= JacobianThP(x);


M=((F*C0*F'));
ML=chol(M,'lower');
Minv=(ML^-1)'*(ML^-1);
C=C0-(C0*F'*(Minv)*F*C0);

% extract just the errors in x and x0

xerr=sqrt(diag(C));
xerr=xerr.*x0mean;
x0err=sqrt(diag(C0));
x0err=x0err.*x0mean;

% actual values of elements x0 and x, and posterior error covariance matrix
x0t=x0.*x0mean;
x=x.*x0mean;
Covx=zeros(length(x),length(x));
for i=1:length(x);
    for j=1:length(x);
        Covx(i,j)=C(i,j)*x0mean(i)*x0mean(j);
    end
end
%%

% plot results

[xresids,fresids]=ThPplots(x,xerr,x0t,x0err,x0mean);

save('stationresults.mat','x','x0t','xerr','x0err','Covx','Thdq','xresids');
