function f = pondcstrCostFcn(X,U,e,data)
% custom cost function 
p = data.PredictionHorizon;
U1 = U(1:p,data.MVIndex(1));
X1 = X(2:p+1,1);
X2 = X(2:p+1,2);

f = sum((X2.*U1.*sqrt(2*9.81*X1)).^2)+ sum((U1.*sqrt(2*9.81*X1)-1/p*sum(U1.*sqrt(2*9.81*X1))).^2);
%f = sum((X2.*U1.*sqrt(2*9.81*X1)).^2);
%f = sum((X2.*U1.*sqrt(2*9.81*X1)).^2)+ sum(U1.*sqrt(2*9.81*X1)-1/p*sum(U1.*sqrt(2*9.81*X1))) + 100*X(end,1);