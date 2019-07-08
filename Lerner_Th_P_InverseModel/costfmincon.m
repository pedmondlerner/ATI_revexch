% This is the cost function for fmincon
% fun is the value of the cost function, and g is its gradient

 function [fun,g]=costfmincon(x);
 
 global x0 C0 

     
 fun=(x-x0)'*inv(C0)*(x-x0);
 
 g=2*inv(C0)*(x-x0);
 
