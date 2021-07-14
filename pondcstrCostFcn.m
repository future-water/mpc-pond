function f = pondcstrCostFcn(X,U,e,data)

%a custom cost function 
% U1 = U(1:end-1,1);
% X1 = X(2:end,1);
% X2 = X(2:end,2);
% 
% f = -data.Ts*sum(sum(X2.*U1.*sqrt(X1)));
% U1 = U(end-1,1);
% X1 = X(end,1);
% X2 = X(end,2); 
% f = -X2.*U1.*sqrt(X1);

% p = data.PredictionHorizon;
% U1 = U(1:p,data.MVIndex(1));
% X1 = X(2:p+1,1);
% X2 = X(2:p+1,2);
% f = sum(sum(X2.*U1.*sqrt(X1)));

p = data.PredictionHorizon;
U1 = U(1:p,data.MVIndex(1));
X1 = X(2:p+1,1);
X2 = X(2:p+1,2);
f = sum(X2.*U1.*sqrt(2*9.81*X1));