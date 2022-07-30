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

f = sum((X2.*q_out).^2) + 10*sum((q_out-1/p*sum(q_out)).^2)+ 900*(X1(end).^2);  % costfun_2
%f = sum((X2.*q_out).^2) + 10*sum((q_out-1/p*sum(q_out)).^2)+ 10*(sum(U1.^2)) + 1000*(X1(end).^2);
%f = sum((X2.*q_out).^2) + 5*sum((q_out-1/p*sum(q_out)).^2) + 10*(X1(end).^2);


%f = sum((X2.*X3).^2) + 1000*sum((X3-1/p*sum(X3)).^2)+10*sum(X1.^2);
%f = sum((X2.*q_out).^2) + 1000*sum((diff(q_out)).^2)+ 10*sum(X1.^2);