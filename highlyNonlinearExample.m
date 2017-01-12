%Ejemplo complejo no lineal
%Example of page 285
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
%--------------Simulation--------------%

%Parameters for the simulation
duration = 100; %seconds
dt = 1; %seconds
n_iter = floor( duration/dt );
t = dt*(0:(n_iter-1));

%Process model
Rww = 10; %Process noise
Rvv = 1; %Measurement noise

%Simulation Arrays
x = zeros(n_iter,1); %true state
y = zeros(n_iter,1); %measurement
%y_true = zeros(n_iter,1); %true measurement

%Real initial condition
x(1) = 0.1; %state
y(1) = (x(1)^2)/20 + sqrt(Rvv)*randn;
%y_true(1) = (x(1)^2)/20;
%Simulation
for k = 2:n_iter
    x(k) = x(k-1)/2 + 25*x(k-1)/(1+x(k-1)^2) + 8*cos(1.2*(k-2)) + sqrt(Rww)*randn;
    %y_true(k) = (x(k)^2)/20;
    y(k) = (x(k)^2)/20 + sqrt(Rvv)*randn;
end

%Show the results of the simulation
figure
subplot(2,1,1)
plot(t,x);
grid, xlabel('Time [s]'), ylabel('State'), title('Simulated State');
subplot(2,1,2)
plot(t,y);
%hold on
%plot(t,y_true);
grid, xlabel('Time [s]'), ylabel('Measurement'), title('Simulated Measurements');
%legend('Measurement','True value');


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
x_est(1) = 0.1; %original 0.1
P_est(1) = 5; %original 5

%Noise
Rww_fil = Rww;
Rvv_fil = Rvv ;

y_pred(1) = (x_est(1)^2)/20;
for k = 2:n_iter
    
    %Prediction
    x_pred(k) = x_est(k-1)/2 + 25*x_est(k-1)/(1+x_est(k-1)^2) + 8*cos(1.2*(k-2));
    y_pred(k) = (x_pred(k)^2)/20;
    A = 0.5 + 25*(1-x_pred(k)^2)/(1+x_pred(k)^2);
    C = x_pred(k)/10;
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
subplot(2,1,1)
plot(t,x);
hold on
plot(t,x_est);
ylim([-30,30])
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('State Estimation with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
legend('True State','EKF Estimate');


subplot(2,1,2)
error = x-x_est;
plot(t,error);
hold on
plot(t,2*sqrt(P_est));
plot(t,-2*sqrt(P_est));
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('State Estimation Error');
legend('Error', '$$+2\sqrt{P_{est}} $$','$$-2\sqrt{P_{est}} $$');
            
%%

%--------------PF Bootstrap---------------
n_part = 100;

%Noise
Rww_fil = Rww;
Rvv_fil = Rvv;

%initial condition
x0 = 0.1;
P0 = 5;

%Inicializacion de particulas
particle = x0 + sqrt(P0)*randn(n_part,1);
particle_pred = zeros(n_part,1);
weight = ones(n_part,1)/n_part;

%Estimadores
x_est_bpf = zeros(n_iter,1);
x_est_bpf(1) = mean(particle);

for k=2:n_iter

    for i=1:n_part
        
        %Importance sampling (prediccion desde k-1 hacia k)
        particle_pred(i) = particle(i)/2 + 25*particle(i)/(1+particle(i)^2) + 8*cos(1.2*(k-2)) + sqrt(Rww_fil)*randn;
        
        %Weight update (medicion de valor en k)
        innov = y(k) - (particle_pred(i)^2)/20;
        weight(i) = exp( -log(sqrt(2*pi*Rvv_fil)) -(( innov )^2)/(2*Rvv_fil) );
    end
    disp(sum(weight));
    %Weigth normalization
    weight = weight/sum(weight);

    %Resampling
    cdf = cumsum(weight);
    for i=1:n_part
        sam = rand;
        ind = find( sam<=cdf ,1);
        particle(i) = particle_pred(ind);
    end
    
    %Estimacion del estado
    x_est_bpf(k) = mean(particle);
end

figure
subplot(2,1,1)
plot(t,x);
hold on
plot(t,x_est_bpf);
ylim([-30,30])
grid, xlabel('Time [s]'), ylabel('State'), title('State Estimation with $$x_0\sim \mathcal{N}(0.1,\ 5)$$');
legend('True State','BPF Estimate');


subplot(2,1,2)
error = x-x_est_bpf;
plot(t,error);
%hold on
%plot(t,2*sqrt(P_est));
%plot(t,-2*sqrt(P_est));
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('State Estimation Error');
%legend('Error', '$$+2\sqrt{P_{est}} $$','$$-2\sqrt{P_{est}} $$');

%%

%Comparison BPF and EKF
figure
plot(t,x,'LineWidth',2,'Color','k');
hold on
plot(t,x_est,'LineWidth',1,'Color','b')
plot(t,x_est_bpf,'LineWidth',1,'Color','r');
ylim([-30,30])
grid, xlabel('Time [s]'), ylabel('State'), title('State Estimation Comparison between EKF and BPF');
legend('True State','EKF Estimate','BPF Estimate');


%%

%RMSE comparison
n_simulation = 30;

rmse_ekf = zeros(n_simulation,1);
rmse_bpf = zeros(n_simulation,1);

for n = 1:n_simulation
    
    %Parameters for the simulation
    duration = 100; %seconds
    dt = 1; %seconds
    n_iter = floor( duration/dt );
    t = dt*(0:(n_iter-1));

    %Process model
    Rww = 10; %Process noise
    Rvv = 1; %Measurement noise

    %Simulation Arrays
    x = zeros(n_iter,1); %true state
    y = zeros(n_iter,1); %measurement
    %y_true = zeros(n_iter,1); %true measurement

    %Real initial condition
    x(1) = 0.1; %state
    y(1) = (x(1)^2)/20 + sqrt(Rvv)*randn;
    %y_true(1) = (x(1)^2)/20;
    %Simulation
    for k = 2:n_iter
        x(k) = x(k-1)/2 + 25*x(k-1)/(1+x(k-1)^2) + 8*cos(1.2*(k-2)) + sqrt(Rww)*randn;
        %y_true(k) = (x(k)^2)/20;
        y(k) = (x(k)^2)/20 + sqrt(Rvv)*randn;
    end
    
    %EKF ESTIMATION
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
    x_est(1) = 0.1; %original 0.1
    P_est(1) = 5; %original 5

    %Noise
    Rww_fil = Rww;
    Rvv_fil = Rvv ;

    y_pred(1) = (x_est(1)^2)/20;
    for k = 2:n_iter

        %Prediction
        x_pred(k) = x_est(k-1)/2 + 25*x_est(k-1)/(1+x_est(k-1)^2) + 8*cos(1.2*(k-2));
        y_pred(k) = (x_pred(k)^2)/20;
        A = 0.5 + 25*(1-x_pred(k)^2)/(1+x_pred(k)^2);
        C = x_pred(k)/10;
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
    
    error_ekf = x-x_est;
    rmse_ekf(n) = sqrt(sum(((x-x_est).^2))/n_iter);
    
    
    %BPF ESTIMATION
    n_part = 100;

    %Noise
    Rww_fil = Rww;
    Rvv_fil = Rvv;

    %initial condition
    x0 = 0.1;
    P0 = 5;

    %Inicializacion de particulas
    particle = x0 + sqrt(P0)*randn(n_part,1);
    particle_pred = zeros(n_part,1);
    weight = ones(n_part,1)/n_part;

    %Estimadores
    x_est_bpf = zeros(n_iter,1);
    x_est_bpf(1) = mean(particle);

    for k=2:n_iter

        for i=1:n_part

            %Importance sampling (prediccion desde k-1 hacia k)
            particle_pred(i) = particle(i)/2 + 25*particle(i)/(1+particle(i)^2) + 8*cos(1.2*(k-2)) + sqrt(Rww_fil)*randn;

            %Weight update (medicion de valor en k)
            innov = y(k) - (particle_pred(i)^2)/20;
            weight(i) = exp( -log(sqrt(2*pi*Rvv_fil)) -(( innov )^2)/(2*Rvv_fil) );
        end

        %Weigth normalization
        weight = weight/sum(weight);

        %Resampling
        cdf = cumsum(weight);
        for i=1:n_part
            sam = rand;
            ind = find( sam<=cdf ,1);
            particle(i) = particle_pred(ind);
        end

        %Estimacion del estado
        x_est_bpf(k) = mean(particle);
    end
    
    
    error_bpf = x-x_est_bpf;
    rmse_bpf(n) = sqrt(sum(((x-x_est_bpf).^2))/n_iter);
    disp('hola');
end

mean_rmse_ekf = mean(rmse_ekf);
std_rmse_ekf = std(rmse_ekf);
mean_rmse_bpf = mean(rmse_bpf);
std_rmse_bpf = std(rmse_bpf);

fprintf('mean rmse ekf: %f \n',mean_rmse_ekf);
fprintf('std rmse ekf: %f \n',std_rmse_ekf); 
fprintf('mean rmse bpf: %f \n',mean_rmse_bpf);
fprintf('std rmse bpf: %f \n',std_rmse_bpf);