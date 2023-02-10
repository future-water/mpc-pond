function f = pondcstrCostFcn(X,U,e,data)
% custom cost function 
p = data.PredictionHorizon;
U1 = U(1:p,data.MVIndex(1));
X1 = X(2:p+1,1);
X2 = X(2:p+1,2);

co = 0.65;
Ao = 1;
g = 32.2;

q_out = U1.*co.*Ao.*sqrt(2*g.*X1).*min(1,X1);

f = 2*sum((X2.*q_out).^2) + 4*sum((q_out-1/p*sum(q_out)).^2)+ 900*(X1(end).^2);
