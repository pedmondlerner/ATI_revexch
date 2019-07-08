% [y]=Jacobian(x);
% This function calculates the jacobian of the model equations (partial first
% derivative of equation in "Thpmodel.m" w.r.t. elements in x). This is used to
% estimate the posterior error covariance matrix.

function y= JacobianThP(x);

% Declare equation normalizing factors as global
global sumTh234d sumTh234p sumTh230d sumTh230p sumTh228d sumTh228p sumP Ddepth Pdepth Thdq cmbTh Th234davg Th234pavg U238avg Th230davg Th230pavg U234avg Th228davg Th228pavg Ra228avg Pavg dTh234 dTh230 dTh228 wavg Bremavg k1avg kdesavg Pdo
dz=diff(Thdq);

% depending on the combination of Th isotopes, define the variables in x
if cmbTh==1;
U238=(x(1:Ddepth));
Th234d=(x(Ddepth+1:2*Ddepth));
Th234p=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
P=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));
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
U234=(x(1:Ddepth));
Th230d=(x(Ddepth+1:2*Ddepth));
Th230p=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1;
P=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

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
Ra228=(x(1:Ddepth));
Th228d=(x(Ddepth+1:2*Ddepth));
Th228p=(x((2*Ddepth)+1:(2*Ddepth)+Pdepth));
if Pdo==1'
P=(x((2*Ddepth)+(Pdepth)+1:(2*Ddepth)+(2*Pdepth)));

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
U238=(x(1:Ddepth));
U234=(x(Ddepth+1:2*Ddepth));
Th234d=(x((2*Ddepth)+1:3*Ddepth));
Th234p=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th230d=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th230p=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

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
U238=(x(1:Ddepth));
Ra228=(x(Ddepth+1:2*Ddepth));
Th234d=(x((2*Ddepth)+1:3*Ddepth));
Th234p=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228d=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228p=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

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
U234=(x(1:Ddepth));
Ra228=(x(Ddepth+1:2*Ddepth));
Th230d=(x((2*Ddepth)+1:3*Ddepth));
Th230p=(x((3*Ddepth)+1:(3*Ddepth)+Pdepth));
Th228d=(x((3*Ddepth)+Pdepth+1:(4*Ddepth)+Pdepth));
Th228p=(x((4*Ddepth)+(Pdepth)+1:(4*Ddepth)+(2*Pdepth)));
if Pdo==1;
P=(x((4*Ddepth)+(2*Pdepth)+1:(4*Ddepth)+(3*Pdepth)));

k1inv=(x((4*Ddepth)+(3*Pdepth)+1:(5*Ddepth)+(3*Pdepth)));
kdesinv=(x((5*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
Breminv=(x((6*Ddepth)+(3*Pdepth)+1:(7*Ddepth)+(3*Pdepth)));
winv=(x((7*Ddepth)+(3*Pdepth)+1:(8*Ddepth)+(3*Pdepth)));
else

k1=(x((4*Ddepth)+(2*Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
kdes=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Brem=(x((6*Ddepth)+(2*Pdepth)+1:(7*Ddepth)+(2*Pdepth)));
w=(x((7*Ddepth)+(2*Pdepth)+1:(8*Ddepth)+(2*Pdepth)));   
end
end    
    
if cmbTh==7;
U238=(x(1:Ddepth));
U234=(x(Ddepth+1:2*Ddepth));
Ra228=(x((2*Ddepth)+1:3*Ddepth));
Th234d=(x((3*Ddepth)+1:(4*Ddepth)));
Th234p=(x((4*Ddepth)+1:(4*Ddepth)+Pdepth));
Th230d=(x((4*Ddepth)+Pdepth+1:(5*Ddepth)+Pdepth));
Th230p=(x((5*Ddepth)+(Pdepth)+1:(5*Ddepth)+(2*Pdepth)));
Th228d=(x((5*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(2*Pdepth)));
Th228p=(x((6*Ddepth)+(2*Pdepth)+1:(6*Ddepth)+(3*Pdepth)));
if Pdo==1;
P=(x((6*Ddepth)+(3*Pdepth)+1:(6*Ddepth)+(4*Pdepth)));

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

% Define the Jacobian depending on the number of Th isotopes used in the
% model.


if Pdo==1;
if cmbTh==1

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(2*Ddepth)+(2*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(3*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(2*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+(Ddepth),j+(3*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+(Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+(Ddepth),j+(5*Ddepth)+(2*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end

if cmbTh==2

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh230*U234avg)/sumTh230d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i,j+(2*Ddepth)+(2*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i,j+(3*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(2*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(Ddepth),j+(3*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(Ddepth),j+(5*Ddepth)+(2*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end

if cmbTh==3

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh228*Ra228avg)/sumTh228d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i,j+(2*Ddepth)+(2*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i,j+(3*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(2*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(Ddepth),j+(3*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(Ddepth),j+(5*Ddepth)+(2*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(2*Ddepth),j+(2*Ddepth)+(Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end


if cmbTh==4

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(3*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(5*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(3*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(5*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(6*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(3*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh230*U234avg)/sumTh230d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(3*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(3*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(4*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(4*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end

if cmbTh==5

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(3*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(5*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(3*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(5*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(6*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(3*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(3*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(3*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(4*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(4*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end



if cmbTh==6

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh230*U234avg)/sumTh230d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i,j+(4*Ddepth)+(3*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i,j+(5*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i,j+(6*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(3*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+Ddepth,j+(5*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+Ddepth,j+(6*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+Ddepth,j+(7*Ddepth)+(3*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(3*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(3*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(4*Ddepth),j+(4*Ddepth)+(2*Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(4*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(4*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end


if cmbTh==7

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;
df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(3*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(4*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(4*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(7*Ddepth)+(4*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(8*Ddepth)+(4*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(3*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(4*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(4*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(6*Ddepth)+(4*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(4*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(8*Ddepth)+(4*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(9*Ddepth)+(4*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh230*U234avg)/sumTh230d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+Pdepth)=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i+(2*Ddepth),j+1+(5*Ddepth)+Pdepth)=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(4*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i+(2*Ddepth),j+(7*Ddepth)+(4*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i+(2*Ddepth),j+(8*Ddepth)+(4*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(3*Ddepth),j+(4*Ddepth)+Pdepth)=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(3*Ddepth),j+1+(5*Ddepth)+Pdepth)=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+Pdepth)=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(3*Ddepth),j+(6*Ddepth)+(4*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(4*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(3*Ddepth),j+(8*Ddepth)+(4*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(3*Ddepth),j+(9*Ddepth)+(4*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(4*Ddepth),j+(2*Ddepth))=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(4*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(4*Ddepth),j+1+(6*Ddepth)+(2*Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(4*Ddepth),j+(6*Ddepth)+(4*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(4*Ddepth),j+(7*Ddepth)+(4*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(4*Ddepth),j+(8*Ddepth)+(4*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(5*Ddepth),j+(5*Ddepth)+(2*Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(5*Ddepth),j+1+(6*Ddepth)+(2*Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(5*Ddepth),j+(6*Ddepth)+(2*Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(5*Ddepth),j+(6*Ddepth)+(4*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(5*Ddepth),j+(7*Ddepth)+(4*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(5*Ddepth),j+(8*Ddepth)+(4*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(5*Ddepth),j+(9*Ddepth)+(4*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));
df(i+(6*Ddepth),j+(6*Ddepth)+(3*Pdepth)+1)=-((Breminv(j)*Bremavg+(winv(j)*wavg/(sumP(j)*dz(j)))))*Pavg;
df(i+(6*Ddepth),j+(6*Ddepth)+(3*Pdepth))=(winv(j)*wavg*Pavg/(sumP(j)*dz(j)));
df(i+(6*Ddepth),j+(8*Ddepth)+(4*Pdepth))=-P(i+1)*Pavg*Bremavg/sumP(j);
df(i+(6*Ddepth),j+(9*Ddepth)+(4*Pdepth))=-((P(i+1)-P(i))*Pavg*wavg/(sumP(j)*dz(i)));
                end
                end
end
end

% if no particles
else
if cmbTh==1

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(2*Ddepth)+(1*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(3*Ddepth)+(1*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(1*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(1*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+(Ddepth),j+(3*Ddepth)+(1*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+(Ddepth),j+(4*Ddepth)+(1*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+(Ddepth),j+(5*Ddepth)+(1*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));

                end
                end
end
end

if cmbTh==2

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh230*U234avg)/sumTh230d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i,j+(2*Ddepth)+(1*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i,j+(3*Ddepth)+(1*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i,j+(4*Ddepth)+(1*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(1*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(Ddepth),j+(3*Ddepth)+(1*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(Ddepth),j+(4*Ddepth)+(1*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(Ddepth),j+(5*Ddepth)+(1*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));

                end
                end
end
end

if cmbTh==3

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;


df(i,j)=(dTh228*Ra228avg)/sumTh228d(j);
df(i,j+(Ddepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i,j+1+(2*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i,j+(2*Ddepth)+(1*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i,j+(3*Ddepth)+(1*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i,j+(4*Ddepth)+(1*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(Ddepth),j+(Ddepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(Ddepth),j+1+(2*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(Ddepth),j+(2*Ddepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(Ddepth),j+(2*Ddepth)+(1*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(Ddepth),j+(3*Ddepth)+(1*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(Ddepth),j+(4*Ddepth)+(1*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(Ddepth),j+(5*Ddepth)+(1*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));

                end
                end
end
end


if cmbTh==4

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(5*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(2*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(5*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(6*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(2*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh230*U234avg)/sumTh230d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(2*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(2*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));

                end
                end
end
end

if cmbTh==5

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(5*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(2*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(2*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(5*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(6*Ddepth)+(2*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(2*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(2*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(2*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));

                end
                end
end
end



if cmbTh==6

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;

df(i,j)=(dTh230*U234avg)/sumTh230d(j);
df(i,j+(2*Ddepth))=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i,j+1+(3*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i,j+(4*Ddepth)+(2*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i,j+(5*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i,j+(6*Ddepth)+(2*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+Ddepth,j+(2*Ddepth))=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+Ddepth,j+1+(3*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+Ddepth,j+(3*Ddepth))=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+Ddepth,j+(4*Ddepth)+(2*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+Ddepth,j+(5*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+Ddepth,j+(6*Ddepth)+(2*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+Ddepth,j+(7*Ddepth)+(2*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(2*Ddepth),j+(3*Ddepth)+(Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(2*Ddepth),j+1+(4*Ddepth)+(Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+(2*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(2*Ddepth),j+(5*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(2*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(3*Ddepth),j+(3*Ddepth)+(Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(3*Ddepth),j+1+(4*Ddepth)+(Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(3*Ddepth),j+(4*Ddepth)+(Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(3*Ddepth),j+(4*Ddepth)+(2*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(3*Ddepth),j+(6*Ddepth)+(2*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(2*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));

                end
                end
end
end


if cmbTh==7

for i=1:(Ddepth);
            for j=1:(Ddepth);
                if i==j;
df(i,j)=(dTh234*U238avg)/sumTh234d(j);
df(i,j+(3*Ddepth))=-(k1inv(j)*k1avg+dTh234)*Th234davg/sumTh234d(j);
df(i,j+1+(4*Ddepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th234pavg/sumTh234d(j);
df(i,j+(6*Ddepth)+(3*Pdepth))=-Th234d(j)*Th234davg*k1avg/sumTh234d(j);
df(i,j+(7*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*kdesavg/sumTh234d(j);
df(i,j+(8*Ddepth)+(3*Pdepth))=Th234p(j+1)*Th234pavg*Bremavg/sumTh234d(j);
df(i+Ddepth,j+(3*Ddepth))=k1inv(j)*k1avg*Th234davg/sumTh234p(j);
df(i+Ddepth,j+1+(4*Ddepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh234+(winv(j)*wavg/dz(j)))*Th234pavg/sumTh234p(j);
df(i+Ddepth,j+(4*Ddepth))=winv(j)*Th234pavg*wavg/(sumTh234p(j)*dz(j));
df(i+Ddepth,j+(6*Ddepth)+(3*Pdepth))=Th234d(j)*Th234davg*k1avg/sumTh234p(j);
df(i+Ddepth,j+(7*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*kdesavg/sumTh234p(j);
df(i+Ddepth,j+(8*Ddepth)+(3*Pdepth))=-Th234p(j+1)*Th234pavg*Bremavg/sumTh234p(j);
df(i+Ddepth,j+(9*Ddepth)+(3*Pdepth))=-(Th234p(j+1)-Th234p(j))*Th234pavg*wavg/(dz(j)*sumTh234p(j));
df(i+(2*Ddepth),j+Ddepth)=(dTh230*U234avg)/sumTh230d(j);
df(i+(2*Ddepth),j+(4*Ddepth)+Pdepth)=-(k1inv(j)*k1avg+dTh230)*Th230davg/sumTh230d(j);
df(i+(2*Ddepth),j+1+(5*Ddepth)+Pdepth)=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th230pavg/sumTh230d(j);
df(i+(2*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-Th230d(j)*Th230davg*k1avg/sumTh230d(j);
df(i+(2*Ddepth),j+(7*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*kdesavg/sumTh230d(j);
df(i+(2*Ddepth),j+(8*Ddepth)+(3*Pdepth))=Th230p(j+1)*Th230pavg*Bremavg/sumTh230d(j);
df(i+(3*Ddepth),j+(4*Ddepth)+Pdepth)=k1inv(j)*k1avg*Th230davg/sumTh230p(j);
df(i+(3*Ddepth),j+1+(5*Ddepth)+Pdepth)=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh230+(winv(j)*wavg/dz(j)))*Th230pavg/sumTh230p(j);
df(i+(3*Ddepth),j+(5*Ddepth)+Pdepth)=winv(j)*Th230pavg*wavg/(sumTh230p(j)*dz(j));
df(i+(3*Ddepth),j+(6*Ddepth)+(3*Pdepth))=Th230d(j)*Th230davg*k1avg/sumTh230p(j);
df(i+(3*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*kdesavg/sumTh230p(j);
df(i+(3*Ddepth),j+(8*Ddepth)+(3*Pdepth))=-Th230p(j+1)*Th230pavg*Bremavg/sumTh230p(j);
df(i+(3*Ddepth),j+(9*Ddepth)+(3*Pdepth))=-(Th230p(j+1)-Th230p(j))*Th230pavg*wavg/(dz(j)*sumTh230p(j));
df(i+(4*Ddepth),j+(2*Ddepth))=(dTh228*Ra228avg)/sumTh228d(j);
df(i+(4*Ddepth),j+(5*Ddepth)+(2*Pdepth))=-(k1inv(j)*k1avg+dTh228)*Th228davg/sumTh228d(j);
df(i+(4*Ddepth),j+1+(6*Ddepth)+(2*Pdepth))=(kdesinv(j)*kdesavg+Breminv(j)*Bremavg)*Th228pavg/sumTh228d(j);
df(i+(4*Ddepth),j+(6*Ddepth)+(3*Pdepth))=-Th228d(j)*Th228davg*k1avg/sumTh228d(j);
df(i+(4*Ddepth),j+(7*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*kdesavg/sumTh228d(j);
df(i+(4*Ddepth),j+(8*Ddepth)+(3*Pdepth))=Th228p(j+1)*Th228pavg*Bremavg/sumTh228d(j);
df(i+(5*Ddepth),j+(5*Ddepth)+(2*Pdepth))=k1inv(j)*k1avg*Th228davg/sumTh228p(j);
df(i+(5*Ddepth),j+1+(6*Ddepth)+(2*Pdepth))=-(kdesinv(j)*kdesavg+Breminv(j)*Bremavg+dTh228+(winv(j)*wavg/dz(j)))*Th228pavg/sumTh228p(j);
df(i+(5*Ddepth),j+(6*Ddepth)+(2*Pdepth))=winv(j)*Th228pavg*wavg/(sumTh228p(j)*dz(j));
df(i+(5*Ddepth),j+(6*Ddepth)+(3*Pdepth))=Th228d(j)*Th228davg*k1avg/sumTh228p(j);
df(i+(5*Ddepth),j+(7*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*kdesavg/sumTh228p(j);
df(i+(5*Ddepth),j+(8*Ddepth)+(3*Pdepth))=-Th228p(j+1)*Th228pavg*Bremavg/sumTh228p(j);
df(i+(5*Ddepth),j+(9*Ddepth)+(3*Pdepth))=-(Th228p(j+1)-Th228p(j))*Th228pavg*wavg/(dz(j)*sumTh228p(j));

                end
                end
end
end
end
    

y=df;