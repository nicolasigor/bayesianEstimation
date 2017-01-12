%Ejemplo sencillo del filtro de Kalman
%Example of page 157
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
%--------------Simulation--------------%

%Parameters for the simulation
duration = 20; %seconds
dt = 100e-3; %seconds
n_iter = floor( duration/dt );
u = ones(n_iter,1)*300e-6; %input
t = dt*(0:(n_iter-1));

%Process model
A = 0.97;
B = 100;
Rww =1e-4; %Process noise

%Measurement model
C = 2;
Rvv = 4; %Measurement noise

%Simulation Arrays
x = zeros(n_iter,1); %true state
y = zeros(n_iter,1); %measurement
y_true = zeros(n_iter,1); %true measurement

%Real initial condition
x(1) = 2.5; %state
y(1) = C*x(1) + sqrt(Rvv)*randn;
y_true(1) = C*x(1);
%Simulation
for k = 2:n_iter
    x(k) = A*x(k-1) + B*u(k-1) + sqrt(Rww)*randn;
    y_true(k) = C*x(k);
    y(k) = C*x(k) + sqrt(Rvv)*randn;
end

%Show the results of the simulation
figure
subplot(2,1,1)
plot(t,x);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('Simulated State with $$R_{ww}=$$ %2.1d',Rww));
subplot(2,1,2)
plot(t,y);
hold on
plot(t,y_true);
grid, xlabel('Time [s]'), ylabel('Measurement'), title(sprintf('Simulated Measurements with $$R_{vv}=$$ %2.1d',Rvv));
legend('Measurement','True value');


%%

%---------------Filtering----------------%

%Filter Arrays
x_est   = zeros(n_iter,1); %state estimate
x_pred  = zeros(n_iter,1); %state estimate prediction
y_pred  = zeros(n_iter,1); %measurement prediction
P_est   = zeros(n_iter,1); %state covariance estimate
P_pred  = zeros(n_iter,1); %state covariance prediction
innov   = zeros(n_iter,1); %innovation
R_innov = zeros(n_iter,1); %innovation covariance
Kgain   = zeros(n_iter,1); %gain
winSm   = zeros(n_iter,1); %average window


%Filter initial conditions
x_est(1) = 0; %original 2.5
P_est(1) = 1; %original 1e-4

%Noise
Rww_fil = Rww;
Rvv_fil = Rvv ;

%size of average window
win = 5; % 2*win es ancho de ventana
winSm(1) = mean(y(1:1+win))/2;

y_pred(1) = C*x_est(1);
for k = 2:n_iter
    
    %Prediction
    x_pred(k) = A*x_est(k-1) + B*u(k-1);
    y_pred(k) = C*x_pred(k);
    P_pred(k) = A*P_est(k-1)*A + Rww_fil;
    
    %Innovation
    innov(k) = y(k) - y_pred(k);
    R_innov(k) = C*P_pred(k)*C + Rvv_fil;
    
    %Kalman Gain
    Kgain(k) = P_pred(k)*C/R_innov(k);
    
    %Update
    x_est(k) = x_pred(k) + Kgain(k)*innov(k);
    P_est(k) = ( 1 - Kgain(k)*C )*P_pred(k);
    
    %Average window
    if k <= win
        winSm(k) = mean(y(1:k+win))/2;
    elseif k < (n_iter-win)
        winSm(k) = mean(y(k-win:k+win))/2;
    else
        winSm(k) = mean(y(k-win:end))/2;
    end
    
end

figure
subplot(2,2,1)
plot(t,x);
hold on
plot(t,x_est);
%plot(t,winSm);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('State Estimation with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
%legend('True State','KF Estimate', 'Average Window');
legend('True State','KF Estimate');


subplot(2,2,2)
error = x-x_est;
plot(t,error);
hold on
plot(t,2*sqrt(P_est));
plot(t,-2*sqrt(P_est));
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('State Estimation Error');
legend('Error', '$$+2\sqrt{P_{est}} $$','$$-2\sqrt{P_{est}} $$');
            
subplot(2,2,3)
plot(t,error);
hold on
plot(t,2*sqrt(P_est));
plot(t,-2*sqrt(P_est));
maximo = max(error);
minimo = min(error);
delta = maximo-minimo;
ylim([minimo-0.25*delta,maximo+0.25*delta]);
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('State Estimation Error');
legend('Error', '$$+2\sqrt{P_{est}} $$','$$-2\sqrt{P_{est}} $$');

subplot(2,2,4)
plot(t(2:end),Kgain(2:end));
hold on
grid, xlabel('Time [s]'), ylabel('Gain');
%title(sprintf('Kalman Gain with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('Kalman Gain Evolution');
%{
figure,
plot(t,x);
hold on
plot(t,x_est);
plot(t,winSm);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('State Estimation with $$x_0=$$ %2.1d and $$P_0=$$ %2.1d',x_est(1), P_est(1)));
legend('True State','KF Estimate', 'Average Window');
%}
