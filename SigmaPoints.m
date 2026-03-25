function [Xi, W] = SigmaPoints(xm,P,kappa)
% kappa: arbitrary constant(normally "n + kappa = 3")
n = numel(xm);        % number of elements: n / xm: X_mean
Xi = zeros(n,2*n+1);  % sigma points = column of Xi
W = zeros(n,1);       % weights

Xi(:,1) = xm;                % Xi_1
W(1) = kappa / (n + kappa);  % W_1
U = chol((n+kappa)*P);       % U'*U = (n+kappa)*P

for k=1:n
    Xi(:,k+1) = xm + U(k,:)'; % u_i = row of U
    W(k+1) = 1/(2*(n+kappa));
end

for k=1:n
    Xi(:,n+k+1) = xm - U(k,:)';
    W(n+k+1) = 1/(2*(n+kappa));
end