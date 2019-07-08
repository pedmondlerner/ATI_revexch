% [c ceq]=Thpmodel(x);
% These are the  thorium and particle dynamic equations.
% When FMINCON finds a solution, the sum of the terms in each equation
% should be 0.
function [c ceq]= Thpmodel(x);

% declare global variables
global sumTh234d sumTh234p sumTh230d sumTh230p sumTh228d sumTh228p sumP Ddepth Pdepth Thdq cmbTh Th234davg Th234pavg U238avg Th230davg Th230pavg U234avg Th228davg Th228pavg Ra228avg Pavg dTh234 dTh230 dTh228 wavg Bremavg k1avg kdesavg Pdo

% vector of depth differences

dzTh=diff(Thdq);
if cmbTh==1;
U238inv=(x(1:Ddepth));
Th234dinv=(x(Ddepth+1:2*Ddepth));
Th234pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));
k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end


if cmbTh==2;
U234inv=(x(1:Ddepth));
Th230dinv=(x(Ddepth+1:2*Ddepth));
Th230pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end

if cmbTh==3;
Ra228inv=(x(1:Ddepth));
Th228dinv=(x(Ddepth+1:2*Ddepth));
Th228pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1'
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
else
k1inv=(x((2*Ddepth)+(Pdepth)+1:(3*Ddepth)+(Pdepth)));
kdesinv=(x((3*Ddepth)+(Pdepth)+1:(4*Ddepth)+(Pdepth)));
Breminv=(x((4*Ddepth)+(Pdepth)+1:(5*Ddepth)+(Pdepth)));
winv=(x((5*Ddepth)+(Pdepth)+1:(6*Ddepth)+(Pdepth)));
end
end
if cmbTh==4;
U238inv=(x(1:Ddepth));
U234inv=(x(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th230dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th230pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 


if cmbTh==5;
U238inv=(x(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else
k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));
end
end 

if cmbTh==6;
U234inv=(x(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));
Th230dinv=(x((2*Ddepth)+1:3*Ddepth));
Th230pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else

k1inv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdesinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Breminv=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
winv=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));   
end
end    
    
if cmbTh==7;
U238inv=(x(1:Ddepth));
U234inv=(x(Ddepth+1:2*Ddepth));
Ra228inv=(x((2*Ddepth)+1:3*Ddepth));
Th234dinv=(x((3*Ddepth)+1:(4*Ddepth)));
Th234pinv=(x((4*Ddepth)+1:(4*Ddepth)+Pdepth));
Th230dinv=(x((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));
Th230pinv=(x((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
Th228dinv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Th228pinv=(x((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
if Pdo==1;
Pinv=(x((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));

k1inv=(x((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));
kdesinv=(x((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));
Breminv=(x((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));
winv=(x((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));
else
k1inv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
kdesinv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
Breminv=(x((8*Ddepth)+(3*Pdepth)+1:(9*Ddepth)+(3*Pdepth)));
winv=(x((9*Ddepth)+(3*Pdepth)+1:(10*Ddepth)+(3*Pdepth)));
end
end

if cmbTh==1
f1=zeros(Ddepth,1);f2=f1;f7=f1;
for j=1:Ddepth
f1(j)=(((U238inv(j)*U238avg*dTh234))-((k1avg*k1inv(j)*Th234dinv(j)*Th234davg))-((dTh234*Th234dinv(j)*Th234davg))+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)))/sumTh234d(j);
f2(j)=((((k1avg*k1inv(j))*Th234dinv(j)*Th234davg))-((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))-((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg))-((dTh234*Th234pinv(j+1)*Th234pavg))-(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))))/sumTh234p(j);     
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
end
f=vertcat(f1,f2,f7);
end

if cmbTh==2
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth     
f3(j)=(((U234inv(j)*U234avg*dTh230))-((k1avg*k1inv(j)*Th230dinv(j)*Th230davg))-((dTh230*Th230dinv(j)*Th230davg))+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)))/sumTh230d(j);
f4(j)=((((k1avg*k1inv(j))*Th230dinv(j)*Th230davg))-((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))-((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg))-((dTh230*Th230pinv(j+1)*Th230pavg))-(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))))/sumTh230p(j);
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
end
f=vertcat(f3,f4,f7);
end

if cmbTh==3
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
f5(j)=(((Ra228avg*Ra228inv(j)*dTh228))-((k1avg*k1inv(j)*Th228dinv(j)*Th228davg))-((dTh228*Th228dinv(j)*Th228davg))+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)))/sumTh228d(j);
f6(j)=((((k1avg*k1inv(j))*Th228dinv(j)*Th228davg))-((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))-((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg))-((dTh228*Th228pinv(j+1)*Th228pavg))-(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))))/sumTh228p(j);
end
f=vertcat(f5,f6,f7);
end

if cmbTh==4
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth
f1(j)=(((U238inv(j)*U238avg*dTh234))-((k1avg*k1inv(j)*Th234dinv(j)*Th234davg))-((dTh234*Th234dinv(j)*Th234davg))+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)))/sumTh234d(j);
f2(j)=((((k1avg*k1inv(j))*Th234dinv(j)*Th234davg))-((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))-((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg))-((dTh234*Th234pinv(j+1)*Th234pavg))-(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))))/sumTh234p(j);      
f3(j)=(((U234inv(j)*U234avg*dTh230))-((k1avg*k1inv(j)*Th230dinv(j)*Th230davg))-((dTh230*Th230dinv(j)*Th230davg))+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)))/sumTh230d(j);
f4(j)=((((k1avg*k1inv(j))*Th230dinv(j)*Th230davg))-((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))-((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg))-((dTh230*Th230pinv(j+1)*Th230pavg))-(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))))/sumTh230p(j);
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
end
f=vertcat(f1,f2,f3,f4,f7);
end
if cmbTh==5
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth
f1(j)=(((U238inv(j)*U238avg*dTh234))-((k1avg*k1inv(j)*Th234dinv(j)*Th234davg))-((dTh234*Th234dinv(j)*Th234davg))+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)))/sumTh234d(j);
f2(j)=((((k1avg*k1inv(j))*Th234dinv(j)*Th234davg))-((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))-((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg))-((dTh234*Th234pinv(j+1)*Th234pavg))-(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))))/sumTh234p(j);      
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
f5(j)=(((Ra228avg*Ra228inv(j)*dTh228))-((k1avg*k1inv(j)*Th228dinv(j)*Th228davg))-((dTh228*Th228dinv(j)*Th228davg))+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)))/sumTh228d(j);
f6(j)=((((k1avg*k1inv(j))*Th228dinv(j)*Th228davg))-((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))-((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg))-((dTh228*Th228pinv(j+1)*Th228pavg))-(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))))/sumTh228p(j);
end
f=vertcat(f1,f2,f5,f6,f7);
end

if cmbTh==6
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth     
f3(j)=(((U234inv(j)*U234avg*dTh230))-((k1avg*k1inv(j)*Th230dinv(j)*Th230davg))-((dTh230*Th230dinv(j)*Th230davg))+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)))/sumTh230d(j);
f4(j)=((((k1avg*k1inv(j))*Th230dinv(j)*Th230davg))-((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))-((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg))-((dTh230*Th230pinv(j+1)*Th230pavg))-(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))))/sumTh230p(j);
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
f5(j)=(((Ra228avg*Ra228inv(j)*dTh228))-((k1avg*k1inv(j)*Th228dinv(j)*Th228davg))-((dTh228*Th228dinv(j)*Th228davg))+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)))/sumTh228d(j);
f6(j)=((((k1avg*k1inv(j))*Th228dinv(j)*Th228davg))-((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))-((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg))-((dTh228*Th228pinv(j+1)*Th228pavg))-(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))))/sumTh228p(j);
end
f=vertcat(f3,f4,f5,f6,f7);
end


if cmbTh==7
f1=zeros(Ddepth,1);f2=f1;f3=f2;f4=f1;f5=f1;f6=f1;f7=f1;
for j=1:Ddepth
f1(j)=(((U238inv(j)*U238avg*dTh234))-((k1avg*k1inv(j)*Th234dinv(j)*Th234davg))-((dTh234*Th234dinv(j)*Th234davg))+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)))/sumTh234d(j);
f2(j)=((((k1avg*k1inv(j))*Th234dinv(j)*Th234davg))-((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg))-((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg))-((dTh234*Th234pinv(j+1)*Th234pavg))-(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))))/sumTh234p(j);      
f3(j)=(((U234inv(j)*U234avg*dTh230))-((k1avg*k1inv(j)*Th230dinv(j)*Th230davg))-((dTh230*Th230dinv(j)*Th230davg))+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)))/sumTh230d(j);
f4(j)=((((k1avg*k1inv(j))*Th230dinv(j)*Th230davg))-((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg))-((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg))-((dTh230*Th230pinv(j+1)*Th230pavg))-(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))))/sumTh230p(j);
if Pdo==1;
f7(j)=(-((Bremavg*Breminv(j)*Pinv(j+1)*Pavg))-(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))))/sumP(j);
else
f7=[];
end
f5(j)=(((Ra228avg*Ra228inv(j)*dTh228))-((k1avg*k1inv(j)*Th228dinv(j)*Th228davg))-((dTh228*Th228dinv(j)*Th228davg))+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)))/sumTh228d(j);
f6(j)=((((k1avg*k1inv(j))*Th228dinv(j)*Th228davg))-((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg))-((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg))-((dTh228*Th228pinv(j+1)*Th228pavg))-(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))))/sumTh228p(j);
end
f=vertcat(f1,f2,f3,f4,f5,f6,f7);
end

ceq=f;

c=[];
