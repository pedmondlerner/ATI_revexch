%%
% [xresids,fresids]=ThPplots(x,xerr,x0,x0err);
% Plots the posterior estimates of the radiochemical activities, particle
% concentration, and rate paramters. Also plots the normalized residuals
% and equations residuals in order to assess, respectively, the fit the of the model to
% the data, and if the model equations are satisfied (i.e. f(x)=0).
%
%
%%%%

function [xresids, fresids]= ThPplots(x,xerr,x0,x0err,x0mean);

% Declare global variables

global cmbTh Pdepth Ddepth Thdq Thdq2 Pdo

% Define the isotope activities, particle concentration, and rate
% paramters, and their errors, for both the prior and posterior estimates
% of x
set(0,'DefaultAxesFontSize',20)


if cmbTh==1;
U238=(x0(1:Ddepth));U238err=(x0err(1:Ddepth));
Th234d=(x0(Ddepth+1:2*Ddepth));Th234derr=(x0err(Ddepth+1:2*Ddepth));
Th234p=(x0((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th234perr=(x0err((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
P=(x0((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Perr=(x0err((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));
k1=(x0((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1err=(x0err((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdes=(x0((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdeserr=(x0err((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Brem=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Bremerr=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
w=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));werr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1=(x0((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1err=(x0err((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdes=(x0((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdeserr=(x0err((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Brem=(x0((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Bremerr=(x0err((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
w=(x0((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));werr=(x0err((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));    
end
end


if cmbTh==2;
U234=(x0(1:Ddepth));U234err=(x0err(1:Ddepth));
Th230d=(x0(Ddepth+1:2*Ddepth));Th230derr=(x0err(Ddepth+1:2*Ddepth));
Th230p=(x0((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th230perr=(x0err((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
P=(x0((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Perr=(x0err((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1=(x0((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1err=(x0err((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdes=(x0((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdeserr=(x0err((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Brem=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Bremerr=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
w=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));werr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1=(x0((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1err=(x0err((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdes=(x0((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdeserr=(x0err((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Brem=(x0((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Bremerr=(x0err((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
w=(x0((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));werr=(x0err((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));    
end
end

if cmbTh==3;
Ra228=(x0(1:Ddepth));Ra228err=(x0err(1:Ddepth));
Th228d=(x0(Ddepth+1:2*Ddepth));Th228derr=(x0err(Ddepth+1:2*Ddepth));
Th228p=(x0((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th228perr=(x0err((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
P=(x0((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Perr=(x0err((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1=(x0((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1err=(x0err((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdes=(x0((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdeserr=(x0err((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Brem=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Bremerr=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
w=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));werr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1=(x0((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1err=(x0err((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdes=(x0((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdeserr=(x0err((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Brem=(x0((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Bremerr=(x0err((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
w=(x0((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));werr=(x0err((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end

if cmbTh==4;
Ra228=(x0(1:Ddepth));Ra228err=(x0err(1:Ddepth));
U234=(x0(Ddepth+1:2*Ddepth));U234err=(x0err(Ddepth+1:2*Ddepth));
Th234d=(x0((2*Ddepth)+1:3*Ddepth));Th234derr=(x0err((2*Ddepth)+1:3*Ddepth));
Th234p=(x0((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th234perr=(x0err((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th230d=(x0((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th230derr=(x0err((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th230p=(x0((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th230perr=(x0err((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x0((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Perr=(x0err((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1=(x0((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1err=(x0err((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdes=(x0((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdeserr=(x0err((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Brem=(x0((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Bremerr=(x0err((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
w=(x0((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));werr=(x0err((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1err=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdes=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdeserr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Brem=(x0((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Bremerr=(x0err((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
w=(x0((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));werr=(x0err((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 


if cmbTh==5;
U238=(x0(1:Ddepth));U238err=(x0err(1:Ddepth));
Ra228=(x0(Ddepth+1:2*Ddepth));Ra228err=(x0err(Ddepth+1:2*Ddepth));
Th234d=(x0((2*Ddepth)+1:3*Ddepth));Th234derr=(x0err((2*Ddepth)+1:3*Ddepth));
Th234p=(x0((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th234perr=(x0err((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228d=(x0((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th228derr=(x0err((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228p=(x0((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th228perr=(x0err((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x0((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Perr=(x0err((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1=(x0((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1err=(x0err((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdes=(x0((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdeserr=(x0err((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Brem=(x0((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Bremerr=(x0err((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
w=(x0((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));werr=(x0err((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1err=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdes=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdeserr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Brem=(x0((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Bremerr=(x0err((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
w=(x0((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));werr=(x0err((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 

if cmbTh==6;
U234=(x0(1:Ddepth));U234err=(x0err(1:Ddepth));
Ra228=(x0(Ddepth+1:2*Ddepth));Ra228err=(x0err(Ddepth+1:2*Ddepth));
Th230d=(x0((2*Ddepth)+1:3*Ddepth));Th230derr=(x0err((2*Ddepth)+1:3*Ddepth));
Th230p=(x0((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th230perr=(x0err((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228d=(x0((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th228derr=(x0err((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228p=(x0((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th228perr=(x0err((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x0((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Perr=(x0err((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1=(x0((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1err=(x0err((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdes=(x0((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdeserr=(x0err((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Brem=(x0((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Bremerr=(x0err((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
w=(x0((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));werr=(x0err((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1=(x0((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1err=(x0err((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdes=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdeserr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Brem=(x0((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Bremerr=(x0err((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
w=(x0((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));werr=(x0err((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end    
    
if cmbTh==7;
U238=(x0(1:Ddepth));U238err=(x0err(1:Ddepth));
U234=(x0(Ddepth+1:2*Ddepth));U234err=(x0err(Ddepth+1:2*Ddepth));
Ra228=(x0((2*Ddepth)+1:3*Ddepth));Ra228err=(x0err((2*Ddepth)+1:3*Ddepth));
Th234d=(x0((3*Ddepth)+1:(4*Ddepth)));Th234derr=(x0err((3*Ddepth)+1:(4*Ddepth)));
Th234p=(x0((4*Ddepth)+1:(4*Ddepth)+Pdepth));Th234perr=(x0err((4*Ddepth)+1:(4*Ddepth)+Pdepth));
Th230d=(x0((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));Th230derr=(x0err((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));
Th230p=(x0((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Th230perr=(x0err((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
Th228d=(x0((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));Th228derr=(x0err((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Th228p=(x0((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));Th228perr=(x0err((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
if Pdo==1;
P=(x0((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));Perr=(x0err((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));

k1=(x0((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));k1err=(x0err((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));
kdes=(x0((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));kdeserr=(x0err((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));
Brem=(x0((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));Bremerr=(x0err((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));
w=(x0((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));werr=(x0err((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));
else
k1=(x0((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));k1err=(x0err((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
kdes=(x0((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));kdeserr=(x0err((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
Brem=(x0((8*Ddepth)+(3*Pdepth)+1:(9*Ddepth)+(3*Pdepth)));Bremerr=(x0err((8*Ddepth)+(3*Pdepth)+1:(9*Ddepth)+(3*Pdepth)));
w=(x0((9*Ddepth)+(3*Pdepth)+1:(10*Ddepth)+(3*Pdepth)));werr=(x0err((9*Ddepth)+(3*Pdepth)+1:(10*Ddepth)+(3*Pdepth)));
end
end


if cmbTh==1;
U238inv=(x(1:Ddepth));U238inverr=(xerr(1:Ddepth));
Th234dinv=(x(Ddepth+1:2*Ddepth));Th234dinverr=(xerr(Ddepth+1:2*Ddepth));
Th234pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th234pinverr=(xerr((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Pinverr=(xerr((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1inverr=(xerr((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Breminverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));winverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1inverr=(xerr((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdesinverr=(xerr((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Breminverr=(xerr((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));winverr=(xerr((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));   
end
end


if cmbTh==2;
U234inv=(x(1:Ddepth));U234inverr=(xerr(1:Ddepth));
Th230dinv=(x(Ddepth+1:2*Ddepth));Th230dinverr=(xerr(Ddepth+1:2*Ddepth));
Th230pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th230pinverr=(xerr((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Pinverr=(xerr((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1inverr=(xerr((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Breminverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));winverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1inverr=(xerr((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdesinverr=(xerr((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Breminverr=(xerr((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));winverr=(xerr((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end

if cmbTh==3;
Ra228inv=(x(1:Ddepth));Ra228inverr=(xerr(1:Ddepth));
Th228dinv=(x(Ddepth+1:2*Ddepth));Th228dinverr=(xerr(Ddepth+1:2*Ddepth));
Th228pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));Th228pinverr=(xerr((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));Pinverr=(xerr((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));k1inverr=(xerr((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Breminverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));winverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));k1inverr=(xerr((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));kdesinverr=(xerr((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));Breminverr=(xerr((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));winverr=(xerr((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end

if cmbTh==4;
Ra228inv=(x(1:Ddepth));Ra228inverr=(xerr(1:Ddepth));
U234inv=(x(Ddepth+1:2*Ddepth));U234inverr=(xerr(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));Th234dinverr=(xerr((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th234pinverr=(xerr((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th230dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th230dinverr=(xerr((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th230pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th230pinverr=(xerr((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Pinverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1inverr=(xerr((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Breminverr=(xerr((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));winverr=(xerr((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1inverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Breminverr=(xerr((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));winverr=(xerr((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 


if cmbTh==5;
U238inv=(x(1:Ddepth));U238inverr=(xerr(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));Ra228inverr=(xerr(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));Th234dinverr=(xerr((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th234pinverr=(xerr((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th228dinverr=(xerr((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th228pinverr=(xerr((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Pinverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1inverr=(xerr((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Breminverr=(xerr((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));winverr=(xerr((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1inverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Breminverr=(xerr((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));winverr=(xerr((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 

if cmbTh==6;
U234inv=(x(1:Ddepth));U234inverr=(xerr(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));Ra228inverr=(xerr(Ddepth+1:2*Ddepth));
Th230dinv=(x((2*Ddepth)+1:3*Ddepth));Th230dinverr=(xerr((2*Ddepth)+1:3*Ddepth));
Th230pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));Th230pinverr=(xerr((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));Th228dinverr=(xerr((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));Th228pinverr=(xerr((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));Pinverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));k1inverr=(xerr((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));Breminverr=(xerr((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));winverr=(xerr((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));k1inverr=(xerr((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));kdesinverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));Breminverr=(xerr((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));winverr=(xerr((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end    
    
if cmbTh==7;
U238inv=(x(1:Ddepth));U238inverr=(xerr(1:Ddepth));
U234inv=(x(Ddepth+1:2*Ddepth));U234inverr=(xerr(Ddepth+1:2*Ddepth));
Ra228inv=(x((2*Ddepth)+1:3*Ddepth));Ra228inverr=(xerr((2*Ddepth)+1:3*Ddepth));
Th234dinv=(x((3*Ddepth)+1:(4*Ddepth)));Th234dinverr=(xerr((3*Ddepth)+1:(4*Ddepth)));
Th234pinv=(x((4*Ddepth)+1:(4*Ddepth)+Pdepth));Th234pinverr=(xerr((4*Ddepth)+1:(4*Ddepth)+Pdepth));
Th230dinv=(x((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));Th230dinverr=(xerr((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));
Th230pinv=(x((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));Th230pinverr=(xerr((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
Th228dinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));Th228dinverr=(xerr((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Th228pinv=(x((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));Th228pinverr=(xerr((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
if Pdo==1;
Pinv=(x((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));Pinverr=(xerr((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));

k1inv=(x((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));k1inverr=(xerr((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));
kdesinv=(x((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));kdesinverr=(xerr((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));
Breminv=(x((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));Breminverr=(xerr((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));
winv=(x((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));winverr=(xerr((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));
else
k1inv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));k1inverr=(xerr((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
kdesinv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));kdesinverr=(xerr((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
Breminv=(x((8*Ddepth)+(3*Pdepth)+1:(9*Ddepth)+(3*Pdepth)));Breminverr=(xerr((8*Ddepth)+(3*Pdepth)+1:(9*Ddepth)+(3*Pdepth)));
winv=(x((9*Ddepth)+(3*Pdepth)+1:(10*Ddepth)+(3*Pdepth)));winverr=(xerr((9*Ddepth)+(3*Pdepth)+1:(10*Ddepth)+(3*Pdepth)));
end
end

% plot activities, concentrations, and rate parameters

if exist('U238inv')==1;
    figure(1);
    errorbar_x(U238inv,Thdq2,U238inverr,U238inverr,'ob');hold on;
    errorbar_x(U238,Thdq2,U238err,U238err, 'xk');xlabel('^{238}U (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('U234inv')==1;
    figure(2);
    errorbar_x(U234inv,Thdq2,U234inverr,U234inverr,'ob');hold on;
    errorbar_x(U234,Thdq2,U234err,U234err,'xk');xlabel('^{234}U (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Ra228inv')==1;
    figure(3);
    errorbar_x(Ra228inv,Thdq2,Ra228inverr,Ra228inverr,'ob');hold on;
    errorbar_x(Ra228,Thdq2,Ra228err,Ra228err,'xk');xlabel('^{228}Ra (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th234dinv')==1;
    figure(4);
    errorbar_x(Th234dinv,Thdq2,Th234dinverr,Th234dinverr,'ob');hold on;
    errorbar_x(Th234d,Thdq2,Th234derr,Th234derr,'xk');xlabel('dissolved ^{234}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th234pinv')==1;
    figure(5);
    errorbar_x(Th234pinv,Thdq,Th234pinverr,Th234pinverr,'ob');hold on;
    errorbar_x(Th234p,Thdq,Th234perr,Th234perr,'xk');xlabel('particulate ^{234}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th230dinv')==1;
    figure(6);
    errorbar_x(Th230dinv,Thdq2,Th230dinverr,Th230dinverr,'ob');hold on;
    errorbar_x(Th230d,Thdq2,Th230derr,Th230derr,'xk');xlabel('dissolved ^{230}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th230pinv')==1;
    figure(7);
    errorbar_x(Th230pinv,Thdq,Th230pinverr,Th230pinverr,'ob');hold on;
    errorbar_x(Th230p,Thdq,Th230perr,Th230perr,'xk');xlabel('particulate ^{230}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th228dinv')==1;
    figure(8);
    errorbar_x(Th228dinv,Thdq2,Th228dinverr,Th228dinverr,'ob');hold on;
    errorbar_x(Th228d,Thdq2,Th228derr,Th228derr,'xk');xlabel('dissolved ^{228}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Th228pinv')==1;
    figure(9);
    errorbar_x(Th228pinv,Thdq,Th228pinverr,Th228pinverr,'ob');hold on;
    errorbar_x(Th228p,Thdq,Th228perr,Th228perr,'xk');xlabel('particulate ^{228}Th (dpm/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end

if exist('Pinv')==1;
    figure(10);
    errorbar_x(Pinv,Thdq,Pinverr,Pinverr,'ob');hold on;
    errorbar_x(P,Thdq,Perr,Perr,'xk');xlabel('Particle Concentration(mg/m^3)');ylabel('depth (m)');set(gca,'ydir','reverse');
end


    figure(11);
    errorbar_x(k1inv,Thdq2,k1inverr,k1inverr,'ob');hold on;
    h=errorbar_x(k1,Thdq2,k1err,k1err,'xk');xlabel('k_1 (yr^{-1})');ylabel('depth (m)');set(gca,'ydir','reverse');set(h,'color',[.6 .6 .6]);
    
    figure(12);
    errorbar_x(kdesinv,Thdq2,kdesinverr,kdesinverr,'ob');hold on;
    h=errorbar_x(kdes,Thdq2,kdeserr,kdeserr,'xk');xlabel('k_{-1} (yr^{-1})');ylabel('depth (m)');set(gca,'ydir','reverse');set(h,'color',[.6 .6 .6]);
    
    figure(13);
    errorbar_x(Breminv,Thdq2,Breminverr,Breminverr,'ob');hold on;
    h=errorbar_x(Brem,Thdq2,Bremerr,Bremerr,'xk');xlabel('\beta_{-1} (yr^{-1})');ylabel('depth (m)');set(gca,'ydir','reverse');set(h,'color',[.6 .6 .6]);
    
    figure(14);
    errorbar_x(winv,Thdq2,winverr,winverr,'ob');hold on;
    h=errorbar_x(w,Thdq2,werr,werr,'xk');xlabel('w (m yr^{-1})');ylabel('depth (m)');set(gca,'ydir','reverse');set(h,'color',[.6 .6 .6]);
    
    % plot normalized residuals. Normalized residuals less than 2 indicate
    % that the model fits the data within +/- standard errors of the data.
    
    xresids=(x-x0)./x0err;
    
    figure(15);
    % the red lines illustrate how many points are within
    % +/-2
    plot(xresids(1:end-(4*Ddepth)),'o');ylabel('normalized residuals');hold on;
    plot(-2.*ones(length(xresids(1:end-(4*Ddepth))),1),'-r');hold on;
    plot(2.*ones(length(xresids(1:end-(4*Ddepth))),1),'-r');
    
    % plot model equation residuals. FMINCON finds a solution,
    % that should satisfy the model equations (f(x)=0) to the 4th order
    % (see constraint tolerance). The largest of these residuals should be no
    % greater than 10^-2.
    
    
    xmodelchk=x./x0mean;
    
    [c f]=Thpmodel(xmodelchk);
    
    if cmbTh==1;
        fnorm=f(:);
    end
    if cmbTh==2;
        fnorm=f(:);
    end
    if cmbTh==3;
        fnorm=f(:);
    end
    if cmbTh==4;
        fnorm=f(:);
    end
    if cmbTh==5;
        fnorm=f(:);
    end
    if cmbTh==6;
        fnorm=f(:);
    end
    if cmbTh==7;
        fnorm=f(:);
    end
    
    figure(16);
    plot(fnorm,'o');ylabel('normalized equation residuals');
    
    fresids=fnorm;
    
    
    