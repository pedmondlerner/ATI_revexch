% This function constructs the prior estimate of the state vector, x0. The
% elements include measured (or interoplated) thorium isotope activities,
% parent activities, particle concentration, and prior estimates of the
% rate constant.


global Th234davg Th234pavg U238avg Th230davg Th230pavg U234avg Th228davg Th228pavg Ra228avg Pavg...
    Th234dnorm Th234pnorm U238norm Th230dnorm Th230pnorm U234norm Th228dnorm Th228pnorm Ra228norm Pnorm...
    Th234derrnorm Th234perrnorm U238errnorm Th230derrnorm Th230perrnorm U234errnorm Th228derrnorm Th228perrnorm Ra228errnorm Perrnorm...
    cmbTh Pdepth Ddepth Thdq Thdq2 k1avg kdesavg Bremavg wavg x0 C0

% these are the prior estimates of the rate constants and their errors
    Brem(1:Ddepth)=1;Bremerr(1:Ddepth)=10;w(1:Ddepth)=700;werr(1:Ddepth)=400;k1(1:Ddepth)=.5;
    k1err(1:Ddepth)=5;kdes(1:Ddepth)=2;kdeserr(1:Ddepth)=5;
    
 k1err=diag(k1err);kdeserr=diag(kdeserr);Bremerr=diag(Bremerr);werr=diag(werr);
 
 % normalize rate constants
 k1avg=mean(k1);kdesavg=mean(kdes);Bremavg=mean(Brem);wavg=mean(w);
 Bremnorm=Brem'./mean(Brem);k1norm=k1'./mean(k1);wnorm=w'./mean(w);kdesnorm=kdes'./mean(kdes);
 Bremerrnorm=(Bremerr./mean(Brem)).^2;k1errnorm=(k1err./mean(k1)).^2;werrnorm=(werr./mean(w)).^2;kdeserrnorm=(kdeserr./mean(kdes)).^2;
 
 % construct vector of means
 if cmbTh==1 || cmbTh==4 || cmbTh==5 || cmbTh==7
 Th234dmean(1:Ddepth,1)=Th234davg;Th234pmean(1:Pdepth,1)=Th234pavg;U238mean(1:Ddepth,1)=U238avg;
 end
 if cmbTh==2 || cmbTh==4 || cmbTh==6 || cmbTh==7
 Th230dmean(1:Ddepth,1)=Th230davg;Th230pmean(1:Pdepth,1)=Th230pavg;U234mean(1:Ddepth,1)=U234avg;
 end
 if cmbTh==3 || cmbTh==5 || cmbTh==6 || cmbTh==7
 Th228dmean(1:Ddepth,1)=Th228davg;Th228pmean(1:Pdepth,1)=Th228pavg;Ra228mean(1:Ddepth,1)=Ra228avg;
 end
 if Pdo==1;
 Pmean(1:Pdepth,1)=Pavg;
 else
 Pnorm=[];Perrnorm=[];Pmean=[];
 end
 k1mean(1:Ddepth,1)=mean(k1);kdesmean(1:Ddepth,1)=mean(kdes);Bremmean(1:Ddepth,1)=mean(Brem);wmean(1:Ddepth,1)=mean(w);
 
 
 
 
 
 % construct prior estimates of state vector and error covariance matrices
 if cmbTh==1
     x0=[U238norm;Th234dnorm;Th234pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U238errnorm,Th234derrnorm,Th234perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U238mean;Th234dmean;Th234pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];
 end
 if cmbTh==2
     x0=[U234norm;Th230dnorm;Th230pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U234errnorm,Th230derrnorm,Th230perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U234mean;Th230dmean;Th230pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];
end
 if cmbTh==3
     x0=[Ra228norm;Th228dnorm;Th228pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(Ra228errnorm,Th228derrnorm,Th228perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[Ra228mean;Th228dmean;Th228pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];

end
 if cmbTh==4
     x0=[U238norm;U234norm;Th234dnorm;Th234pnorm;Th230dnorm;Th230pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U238errnorm,U234errnorm,Th234derrnorm,Th234perrnorm,Th230derrnorm,Th230perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U238mean;U234mean;Th234dmean;Th234pmean;Th230dmean;Th230pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];
 end
 if cmbTh==5
     x0=[U238norm;Ra228norm;Th234dnorm;Th234pnorm;Th228dnorm;Th228pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U238errnorm,Ra228errnorm,Th234derrnorm,Th234perrnorm,Th228derrnorm,Th228perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U238mean;Ra228mean;Th234dmean;Th234pmean;Th228dmean;Th228pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];
 end
 if cmbTh==6
     x0=[U234norm;Ra228norm;Th230dnorm;Th230pnorm;Th228dnorm;Th228pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U234errnorm,Ra228errnorm,Th230derrnorm,Th230perrnorm,Th228derrnorm,Th228perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U234mean;Ra228mean;Th230dmean;Th230pmean;Th228dmean;Th228pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];
 end
 if cmbTh==7
     x0=[U238norm;U234norm;Ra228norm;Th234dnorm;Th234pnorm;Th230dnorm;Th230pnorm;Th228dnorm;Th228pnorm;Pnorm;k1norm;kdesnorm;Bremnorm;wnorm];
     C0=blkdiag(U238errnorm,U234errnorm,Ra228errnorm,Th234derrnorm,Th234perrnorm,Th230derrnorm,Th230perrnorm,Th228derrnorm,Th228perrnorm,Perrnorm,k1errnorm,kdeserrnorm,Bremerrnorm,werrnorm);
     x0mean=[U238mean;U234mean;Ra228mean;Th234dmean;Th234pmean;Th230dmean;Th230pmean;Th228dmean;Th228pmean;Pmean;k1mean;kdesmean;Bremmean;wmean];

 end
