% [sumTh234d,sumTh234p,sumTh230d,sumTh230p,sumTh228d,sumTh228p]=modelnorm(x);
% Outputs the sum of the square of the terms in the thorium and particle
% dynamic equations, to be used as normalization factors in the model.
% This is so that the model equation for each isotope and the
% particles have equal weight in the inversion

function [sumTh234d,sumTh234p,sumTh230d,sumTh230p,sumTh228d,sumTh228p,sumP]=modelcheck(x);

global Th234davg Th234pavg U238avg Th230davg Th230pavg U234avg Th228davg Th228pavg Ra228avg Pavg cmbTh Pdepth Ddepth Thdq Thdq2 dTh234 dTh230 dTh228 wavg Bremavg k1avg kdesavg

% vector of depth differences
dzTh=diff(Thdq);

% depending on the combination of Th isotopes, define the variables in x0
% (note that in mainThp.m, this routine is only run using x0 as an input so
% that the normalization factors are constant).

if cmbTh==1;
U238inv=(x(1:Ddepth));
Th234dinv=(x(Ddepth+1:2*Ddepth));
Th234pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
end


if cmbTh==2;
U234inv=(x(1:Ddepth));
Th230dinv=(x(Ddepth+1:2*Ddepth));
Th230pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
end

if cmbTh==3;
U238inv=(x(1:Ddepth));
Th228dinv=(x(Ddepth+1:2*Ddepth));
Th228pinv=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
Pinv=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

k1inv=(x((2*Ddepth)+(2*Pdepth)+1:(3*Ddepth)+(2*Pdepth)));
kdesinv=(x((3*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Breminv=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
winv=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
end

if cmbTh==4;
Ra228inv=(x(1:Ddepth));
U234inv=(x(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th230dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th230pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
end 


if cmbTh==5;
U238inv=(x(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));
Th234dinv=(x((2*Ddepth)+1:3*Ddepth));
Th234pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
end 

if cmbTh==6;
U234inv=(x(1:Ddepth));
Ra228inv=(x(Ddepth+1:2*Ddepth));
Th230dinv=(x((2*Ddepth)+1:3*Ddepth));
Th230pinv=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228dinv=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228pinv=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
Pinv=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
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
Pinv=(x((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));

k1inv=(x((6*Ddepth)+(4*Pdepth)+1:(7*Ddepth)+(4*Pdepth)));
kdesinv=(x((7*Ddepth)+(4*Pdepth)+1:(8*Ddepth)+(4*Pdepth)));
Breminv=(x((8*Ddepth)+(4*Pdepth)+1:(9*Ddepth)+(4*Pdepth)));
winv=(x((9*Ddepth)+(4*Pdepth)+1:(10*Ddepth)+(4*Pdepth)));
end

% depending on the number of Th isotopes, calculate the sum of the square
% of the terms in teh appropriate model equations.

if cmbTh==1
for j=1:Ddepth
f1(j)=((U238inv(j)*U238avg*dTh234)^2)+((k1avg*k1inv(j)*Th234dinv(j)*Th234davg)^2)+((dTh234*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2);
f2(j)=(((k1avg*k1inv(j))*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2)+((dTh234*Th234pinv(j+1)*Th234pavg)^2)+(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))^2);      
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
end
sumTh234d=sqrt(f1');
sumTh234p=sqrt(f2');
sumTh230d=0;
sumTh230p=0;
sumTh228d=0;
sumTh228p=0;
sumP=sqrt(f7');
end

if cmbTh==2
for j=1:Ddepth      
f3(j)=((U234inv(j)*U234avg*dTh230)^2)+((k1avg*k1inv(j)*Th230dinv(j)*Th230davg)^2)+((dTh230*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2);
f4(j)=(((k1avg*k1inv(j))*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2)+((dTh230*Th230pinv(j+1)*Th230pavg)^2)+(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))^2);
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
end
sumTh234d=0;
sumTh234p=0;
sumTh230d=sqrt(f3');
sumTh230p=sqrt(f4');
sumTh228d=0;
sumTh228p=0;
sumP=sqrt(f7');
end

if cmbTh==3
for j=1:Ddepth
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
f5(j)=((Ra228avg*Ra228inv(j)*dTh228)^2)+((k1avg*k1inv(j)*Th228dinv(j)*Th228davg)^2)+((dTh228*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2);
f6(j)=(((k1avg*k1inv(j))*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2)+((dTh228*Th228pinv(j+1)*Th228pavg)^2)+(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))^2);
end
sumTh234d=0;
sumTh234p=0;
sumTh230d=0;
sumTh230p=0;
sumTh228d=sqrt(f5');
sumTh228p=sqrt(f6');
sumP=sqrt(f7');
end

if cmbTh==4
for j=1:Ddepth
f1(j)=((U238inv(j)*U238avg*dTh234)^2)+((k1avg*k1inv(j)*Th234dinv(j)*Th234davg)^2)+((dTh234*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2);
f2(j)=(((k1avg*k1inv(j))*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2)+((dTh234*Th234pinv(j+1)*Th234pavg)^2)+(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))^2);      
f3(j)=((U234inv(j)*U234avg*dTh230)^2)+((k1avg*k1inv(j)*Th230dinv(j)*Th230davg)^2)+((dTh230*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2);
f4(j)=(((k1avg*k1inv(j))*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2)+((dTh230*Th230pinv(j+1)*Th230pavg)^2)+(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))^2);
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);

end
sumTh234d=sqrt(f1');
sumTh234p=sqrt(f2');
sumTh230d=sqrt(f3');
sumTh230p=sqrt(f4');
sumTh228d=0;
sumTh228p=0;
sumP=sqrt(f7');
end
if cmbTh==5
for j=1:Ddepth
f1(j)=((U238inv(j)*U238avg*dTh234)^2)+((k1avg*k1inv(j)*Th234dinv(j)*Th234davg)^2)+((dTh234*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2);
f2(j)=(((k1avg*k1inv(j))*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2)+((dTh234*Th234pinv(j+1)*Th234pavg)^2)+(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))^2);      
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
f5(j)=((Ra228avg*Ra228inv(j)*dTh228)^2)+((k1avg*k1inv(j)*Th228dinv(j)*Th228davg)^2)+((dTh228*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2);
f6(j)=(((k1avg*k1inv(j))*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2)+((dTh228*Th228pinv(j+1)*Th228pavg)^2)+(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))^2);
end
sumTh234d=sqrt(f1');
sumTh234p=sqrt(f2');
sumTh230d=0;
sumTh230p=0;
sumTh228d=sqrt(f5');
sumTh228p=sqrt(f6');
sumP=sqrt(f7');
end

if cmbTh==6
for j=1:Ddepth     
f3(j)=((U234inv(j)*U234avg*dTh230)^2)+((k1avg*k1inv(j)*Th230dinv(j)*Th230davg)^2)+((dTh230*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2);
f4(j)=(((k1avg*k1inv(j))*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2)+((dTh230*Th230pinv(j+1)*Th230pavg)^2)+(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))^2);
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
f5(j)=((Ra228avg*Ra228inv(j)*dTh228)^2)+((k1avg*k1inv(j)*Th228dinv(j)*Th228davg)^2)+((dTh228*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2);
f6(j)=(((k1avg*k1inv(j))*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2)+((dTh228*Th228pinv(j+1)*Th228pavg)^2)+(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))^2);
end
sumTh234d=0;
sumTh234p=0;
sumTh230d=sqrt(f3');
sumTh230p=sqrt(f4');
sumTh228d=sqrt(f5');
sumTh228p=sqrt(f6');
sumP=sqrt(f7');
end


if cmbTh==7
for j=1:Ddepth
f1(j)=((U238inv(j)*U238avg*dTh234)^2)+((k1avg*k1inv(j)*Th234dinv(j)*Th234davg)^2)+((dTh234*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2);
f2(j)=(((k1avg*k1inv(j))*Th234dinv(j)*Th234davg)^2)+((kdesavg*kdesinv(j)*Th234pinv(j+1)*Th234pavg)^2)+((Bremavg*Breminv(j)*Th234pinv(j+1)*Th234pavg)^2)+((dTh234*Th234pinv(j+1)*Th234pavg)^2)+(((winv(j)*wavg*Th234pavg/dzTh(j))*(Th234pinv(j+1)-Th234pinv(j)))^2);      
f3(j)=((U234inv(j)*U234avg*dTh230)^2)+((k1avg*k1inv(j)*Th230dinv(j)*Th230davg)^2)+((dTh230*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2);
f4(j)=(((k1avg*k1inv(j))*Th230dinv(j)*Th230davg)^2)+((kdesavg*kdesinv(j)*Th230pinv(j+1)*Th230pavg)^2)+((Bremavg*Breminv(j)*Th230pinv(j+1)*Th230pavg)^2)+((dTh230*Th230pinv(j+1)*Th230pavg)^2)+(((winv(j)*wavg*Th230pavg/dzTh(j))*(Th230pinv(j+1)-Th230pinv(j)))^2);
f7(j)=((Bremavg*Breminv(j)*Pinv(j+1)*Pavg)^2)+(((winv(j)*wavg*Pavg/dzTh(j))*((Pinv(j+1))-(Pinv(j))))^2);
f5(j)=((Ra228avg*Ra228inv(j)*dTh228)^2)+((k1avg*k1inv(j)*Th228dinv(j)*Th228davg)^2)+((dTh228*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2);
f6(j)=(((k1avg*k1inv(j))*Th228dinv(j)*Th228davg)^2)+((kdesavg*kdesinv(j)*Th228pinv(j+1)*Th228pavg)^2)+((Bremavg*Breminv(j)*Th228pinv(j+1)*Th228pavg)^2)+((dTh228*Th228pinv(j+1)*Th228pavg)^2)+(((winv(j)*wavg*Th228pavg/dzTh(j))*(Th228pinv(j+1)-Th228pinv(j)))^2);
end
sumTh234d=sqrt(f1');
sumTh234p=sqrt(f2');
sumTh230d=sqrt(f3');
sumTh230p=sqrt(f4');
sumTh228d=sqrt(f5');
sumTh228p=sqrt(f6');
sumP=sqrt(f7');
end