function f = pondcstrCostFcn(X,U,e,data)
% custom cost function 
p = data.PredictionHorizon;
U1 = U(1:p,data.MVIndex(1));
X1 = X(2:p+1,1);
X2 = X(2:p+1,2);

co = 0.67;
Ao = 1.5;
g = 9.81;

qout = U1.*co.*Ao.*sqrt(2*g.*X1);

f = sum((X2.*qout).^2) + 10*sum((qout-1/p*sum(qout)).^2);
