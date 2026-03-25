%% w12 Unscented Kalman Filter Example
%% Attitude estimation by accelerometers & gyroscopes
clear all; close all; clc
d2r = pi/180; r2d = 1/d2r; g = 9.81;
%% Data from gyroscopes & accelerometers
% load real gyroscope data
load ArsGyro
p = wx; q = wy; r = wz;
% load real accelerometer data
load ArsAccel
ax = fx; ay = fy; az = fz;
%% Dynamics
T = 0.01;           % sampling time: taking a picture of a ball every 1 sec
Gamma = eye(3,3);   % Gamma (related to system noise)
H = [
    1 0 0
    0 1 0
    ];              % Measurement matrix: ONLY phi & theta
%% Initialization of the sizes of vectors & matrices
Nsamples = length(p);
t = 0:T:(Nsamples-1)*T;
nx = 3;                    % number of state;
nm = 2;                    % number of measurement
% x = zeros(nx,length(t))
xp = zeros(nx,length(t));
y = zeros(nm,length(t));
sigma_Pp = zeros(nx,length(t));
%% Initial Conditions
% xp = [phi; theta; psi]
xp(:,1) = [0; 0; 0];            % guess of initial posteriori estimation
y(:,1) = H*xp(:,1);             % initialize the measurement (not used)
Pp = 1e3*eye(nx);               % guess of initial error covariance
sigma_Pp(:,1) = sqrt(diag(Pp));
%% Noise
sigma_w = [sqrt(0.0001); sqrt(0.0001); sqrt(0.1)]; % system noise (std of acceleration) ** std = standard deviation(표준편차)
% system noise는 원래 없는게 맞지만 아주 작게 넣었다. 
sigma_v = [sqrt(6); sqrt(6)];                      % measurement noise (std of position sensor)
Q = diag(sigma_w.^2);    % system noise covariance matrix
R = diag(sigma_v.^2);    % measurement noise covariance matrix

%% KF Routine
tic
for i = 1:length(t)-1
    %% True dynamics(실제 데이터를 받아오므로 필요 없음)
    %% Equation 1: Prediction & Equation 2: Prediction of covariance
    %% Sigma points and weights: [Xi,W]=SigmaPoints(xm,P,kappa)
    [Xi, W] = SigmaPoints(xp(:,i), Pp, 0); % kappa = 0 because n+k = 3
    %% ===========================================================
    %% Time update
    %% ===========================================================
    rates =[p(i); q(i); r(i)];
    fXi = zeros(nx, 2*nx+1);
    for k = 1:2*nx+1
        fXi(:,k) = fn_fx(Xi(:,k),rates,T); % 2n+1개의 sigma point를 다음 스텝으로 한 번에 보낸다.
    end
    % Unscented transformation
    [x_, P_] = UT(fXi, W, Q);
    %% ===========================================================
    %% Measurement update
    %% ===========================================================
    %% measurement generation (Equation 3: Innovation Covariance/Equation 4: Residual)
%     Theta = asin(  ax(i)/g );
%     Phi   = asin( -ay(i)/(g*cos(Theta)) );
     Phi = atan(ay(i)/az(i));
     Theta = atan(ax(i)/sqrt(ay(i)^2+az(i)^2));
     y(:,i+1) = [Phi; Theta];
    %% measurement update
     hXi = zeros(nm, 2*nx+1);
     for k = 1:2*nx+1
         hXi(:,k) = fn_hx(fXi(:,k));
     end
     % Unscented transformation
     [yhat, S] = UT(hXi, W, R);
     
     % residual
     nu = y(:,i+1) - yhat;
     
     % for Kalman Gain
     Pxy = zeros(nx,nm);
     for k = 1:2*nx+1
         Pxy = Pxy + W(k)*(fXi(:,k) - x_)*(hXi(:,k) - yhat)';
     end
    %% Equation 5: Kalman gain
    K = Pxy*S^-1;
    %% Equation 6: State update
    xp(:,i+1) = x_ + K*nu;
    %% Equation 7: Covariance update
    Pp = P_ - K*S*K';
    %% ===========================================================
    %% storing error covariance for ploting
    sigma_Pp(:,i+1) = sqrt(diag(Pp));
    %% storing Kalman gain for ploting
    K_store(i) = norm(K);
end
toc
%% Plot: phi
figure(1); clf;
subplot(3,1,1)
plot(t,xp(1,:)*r2d); hold on
plot(t,30*ones(Nsamples,1),':',t,-30*ones(Nsamples,1),':');
xlabel('time, s'); ylabel('\phi, deg')
axis([-inf inf -50 50]); grid
%% Plot: theta
subplot(3,1,2); hold on
plot(t,xp(2,:)*r2d); hold on
plot(t,30*ones(Nsamples,1),':',t,-30*ones(Nsamples,1),':');
xlabel('time, s'); ylabel('\theta, deg')
axis([-inf inf -50 50]); grid
%% Plot: psi
subplot(3,1,3); hold on
plot(t,xp(3,:)*r2d);
xlabel('time, s'); ylabel('\psi, deg')
axis([-inf inf -50 50]); grid

%% Plot: Kalman gain norm
figure(2)
plot(t(1:end-1),K_store); grid on
xlabel('time, s'); ylabel('Kalman gain norm')