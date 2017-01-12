%Ejemplo sencillo no lineal
%Example of page 165, 171
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
%--------------Simulation--------------%

%Parameters for the simulation
duration = 1.5; %seconds
dt = 0.01; %seconds
n_iter = floor( duration/dt );
t = dt*(0:(n_iter-1));

%Process model
Rww = 0; %Process noise
Rvv = 0.09; %Measurement noise

%Simulation Arrays
x = zeros(n_iter,1); %true state
y = zeros(n_iter,1); %measurement
y_true = zeros(n_iter,1); %true measurement

%Real initial condition
x(1) = 2.0; %state
y(1) = x(1)^2 + x(1)^3 + sqrt(Rvv)*randn;
y_true(1) = x(1)^2 + x(1)^3;
%Simulation
for k = 2:n_iter
    x(k) = (1-0.05*dt)*x(k-1) + 0.04*dt*(x(k-1)^2) + sqrt(Rww)*randn;
    y_true(k) = x(k)^2 + x(k)^3;
    y(k) = x(k)^2 + x(k)^3 + sqrt(Rvv)*randn;
end

%Show the results of the simulation
figure
subplot(2,1,1)
plot(t,x);
grid, xlabel('Time [s]'), ylabel('State'), title('Simulated State');
subplot(2,1,2)
plot(t,y);
hold on
plot(t,y_true);
grid, xlabel('Time [s]'), ylabel('Measurement'), title('Simulated Measurements');
legend('Measurement','True value');

%%

%---------------EKF Filtering----------------%

%Filter Arrays
x_est   = zeros(n_iter,1); %state estimate
x_pred  = zeros(n_iter,1); %state estimate prediction
y_pred  = zeros(n_iter,1); %measurement prediction
P_est   = zeros(n_iter,1); %state covariance estimate
P_pred  = zeros(n_iter,1); %state covariance prediction
innov   = zeros(n_iter,1); %innovation
R_innov = zeros(n_iter,1); %innovation covariance
Kgain   = zeros(n_iter,1); %gain

%Filter initial conditions
x_est(1) = 2.0; %original 2.0
P_est(1) = 1e-2; %original 0.01

%Noise
Rww_fil = Rww;
Rvv_fil = Rvv ;

y_pred(1) = x_est(1)^2 + x_est(1)^3;
for k = 2:n_iter
    
    %Prediction
    x_pred(k) = (1-0.05*dt)*x_est(k-1) + 0.04*dt*(x_est(k-1)^2);
    y_pred(k) = x_pred(k)^2 + x_pred(k)^3;
    A = (1-0.05*dt) + 0.08*dt*x_pred(k);
    C = 2*x_pred(k) + 3*(x_pred(k)^2);
    P_pred(k) = A*P_est(k-1)*A + Rww_fil;
    
    %Innovation
    innov(k) = y(k) - y_pred(k);
    R_innov(k) = C*P_pred(k)*C + Rvv_fil;
    
    %Kalman Gain
    Kgain(k) = P_pred(k)*C/R_innov(k);
    
    %Update
    x_est(k) = x_pred(k) + Kgain(k)*innov(k);
    P_est(k) = ( 1 - Kgain(k)*C )*P_pred(k);
        
end

figure
subplot(2,2,1)
plot(t,x);
hold on
plot(t,x_est);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('State Estimation with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
legend('True State','EKF Estimate');


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
ylim([minimo-0*delta,maximo+0*delta]);
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

%%
%---------------EKF Filtering Iterative----------------%

%Filter Arrays
x_est   = zeros(n_iter,4); %state estimate
x_pred  = zeros(n_iter,1); %state estimate prediction
y_pred  = zeros(n_iter,1); %measurement prediction
P_est   = zeros(n_iter,1); %state covariance estimate
P_pred  = zeros(n_iter,1); %state covariance prediction
innov   = zeros(n_iter,1); %innovation
R_innov = zeros(n_iter,1); %innovation covariance
Kgain   = zeros(n_iter,1); %gain

%Filter initial conditions
x_est(1,1) = 2.1; %original 2.0
x_est(1,2) = 2.1; %original 2.0
x_est(1,3) = 2.1; %original 2.0
x_est(1,4) = 2.1; %original 2.0

%Noise
Rww_fil = Rww;
Rvv_fil = Rvv ;

P_initialValues = [1, 1e-2, 1e-4, 1e-6];

for m = 1:4
    P_est(1) = P_initialValues(m);
    y_pred(1) = x_est(1)^2 + x_est(1)^3;
    for k = 2:n_iter

        %Prediction
        x_pred(k) = (1-0.05*dt)*x_est(k-1,m) + 0.04*dt*(x_est(k-1,m)^2);
        y_pred(k) = x_pred(k)^2 + x_pred(k)^3;
        A = (1-0.05*dt) + 0.08*dt*x_pred(k);
        C = 2*x_pred(k) + 3*(x_pred(k)^2);
        P_pred(k) = A*P_est(k-1)*A + Rww_fil;

        %Innovation
        innov(k) = y(k) - y_pred(k);
        R_innov(k) = C*P_pred(k)*C + Rvv_fil;

        %Kalman Gain
        Kgain(k) = P_pred(k)*C/R_innov(k);

        %Update
        x_est(k,m) = x_pred(k) + Kgain(k)*innov(k);
        P_est(k) = ( 1 - Kgain(k)*C )*P_pred(k);
    end
end

figure
plot(t,x);
hold on
plot(t,x_est(:,1));
plot(t,x_est(:,2));
plot(t,x_est(:,3));
plot(t,x_est(:,4));
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('State Estimation with $$x_0=$$ %1.1f',x_est(1,1)));
legend1 = sprintf('EKF Estimate with $$P_0=$$ %1.0d',P_initialValues(1));
legend2 = sprintf('EKF Estimate with $$P_0=$$ %1.0d',P_initialValues(2));
legend3 = sprintf('EKF Estimate with $$P_0=$$ %1.0d',P_initialValues(3));
legend4 = sprintf('EKF Estimate with $$P_0=$$ %1.0d',P_initialValues(4));
legend('True State',legend1, legend2, legend3, legend4);


%%

%--------------PF Bootstrap---------------
n_part = 50;

%Noise
Rww_fil = Rww ;
Rvv_fil = Rvv;

%Resampling
N_t = n_part;

%Filter initial conditions
x0 = 2.1; %original 2.0
P0 = 1e-2; %original 0.01

%Inicializacion de particulas
%particle = x0 + sqrt(P0)*randn(n_part,1);
particle = x0*ones(n_part,1);
particle_pred = zeros(n_part,1);
weight = ones(n_part,1)/n_part;

%Estimadores
x_est_bpf = zeros(n_iter,1);
x_est_bpf(1) = mean(particle);

for k=2:n_iter

    for i=1:n_part
        
        %Importance sampling (prediccion desde k-1 hacia k)
        particle_pred(i) = (1-0.05*dt)*particle(i) + 0.04*dt*(particle(i)^2) + sqrt(Rww_fil)*randn;
        
        %Weight update (medicion de valor en k)
        innov = y(k) - particle_pred(i)^2 + particle_pred(i)^3;
        weight(i) = exp( -log(sqrt(2*pi*Rvv_fil)) -(( innov )^2)/(2*Rvv_fil) );
    end
    %disp(sum(weight));
    %Weigth normalization
    weight = weight/sum(weight);
    
    N_eff = 1/( sum(weight.^2) );
    if N_eff < N_t
        %Resampling
        cdf = cumsum(weight);
        %Systematic resampling
        sam = rand/n_part;
        for i=1:n_part
            samInd = sam + (i-1)/n_part;
            ind = find( samInd<=cdf ,1);
            particle(i) = particle_pred(ind);
        end
    else
        for i=1:n_part
            particle(i) = particle_pred(i);
        end
        
    end
    %Estimacion del estado
    x_est_bpf(k) = mean(particle);
end

figure
subplot(2,1,1)
plot(t,x);
hold on
plot(t,x_est_bpf);
grid, xlabel('Time [s]'), ylabel('State');
title(sprintf('State Estimation with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x0, P0));
legend('True State','BPF Estimate');

subplot(2,1,2)
error = x-x_est_bpf;
plot(t,error);
grid, xlabel('Time [s]'), ylabel('Error');
title('State Estimation Error');