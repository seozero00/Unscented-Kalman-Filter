clear all; clc

xm = 5;     % mean of state x, n = 1
Px = 9;     % covariance of state x
kappa = 2;  % because n+k = 3

[Xi, W] = SigmaPoints(xm,Px,kappa)
[xAvg, xCov] = UT(Xi,W,0)