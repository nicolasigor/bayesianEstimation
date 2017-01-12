% Ejemplo de baterias con mediciones simuladas

set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%Parametros de la bateria (bateria 2 en paper de Pola)
alfa = 0.15;
beta = 12;
gamma = 6.61;
v0 = 4.00;
vl = 3.813;
Ecrit = 19865;
Zp = 0.20;
fprintf('alfa: %f \n',alfa);
fprintf('beta: %f \n',beta);
fprintf('gamma: %f \n',gamma); 
fprintf('v_0: %f \n',v0);
fprintf('v_L: %f \n',vl); 
fprintf('Ecrit: %f \n',Ecrit); 
fprintf('Zp: %f \n',Zp); 

%Leer entrada de corriente
input = csvread('2015_vali1_current.csv');
i_meas = input(:,2);
L = length(i_meas);
dt = 2; %paso de tiempo en segundos
t = dt*(0:(L-1));

%Simular voltaje
v = zeros(L,1);
v_meas = zeros(L,1);
x1 = zeros(L,1);
x2 = zeros(L,1);
%Noise
sigma_1 = 1e-8;
sigma_2 = 1e-8;
sigma_n = 1e-4;
%Initial conditions
x1(1) = Zp;
x2(1) = 1;
v(1) = v0 - i_meas(1)*x1(1);
v_meas(1) = v(1) + sqrt(sigma_n)*randn;

for k = 2:L
    %State transition
    x1(k) = x1(k-1) + sqrt(sigma_1)*randn;
    x2(k) = x2(k-1) - v(k-1)*i_meas(k-1)*dt/Ecrit + sqrt(sigma_2)*randn;
    if x2(k)<0
        x2(k) = 0;
    end
    %Measurement
    v(k) = vl + (v0-vl)*exp(gamma*(x2(k)-1))...
        + alfa*vl*(x2(k)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2(k))))...
        - i_meas(k)*x1(k);
    v_meas(k) =  v(k) + sqrt(sigma_n)*randn;
end

%Mostrar los datos simulados
figure
subplot(2,2,1)
plot(t,i_meas);
grid, xlabel('Time [s]'), ylabel('Current [A]'), title('Measured Discharge Current');
subplot(2,2,2)
plot(t,v_meas);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Simulated Measured Discharge Voltage');

%Mostrar los datos simulados
%figure
subplot(2,2,4)
plot(t,v);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Simulated True Discharge Voltage');
subplot(2,2,3)
plot(t,x2);
grid, xlabel('Time [s]'), ylabel('SOC'), title('Simulated True SOC');

%%
%Filtrado
%--------------PF Bootstrap---------------
n_part = 100;

%Ruidos filtro
w_1 = sigma_1;
w_2 = 1e-5; %sigma_2;
Rnn = 1e-2; %sigma_n;

%Resampling
N_t = n_part;

%Inicializacion de particulas
particle = zeros(n_part,2);
particle_pred =  zeros(n_part,2);
particle(:,1) = ones(n_part,1)*Zp;
particle(:,2) = unifrnd(0.8,0.9, n_part, 1); %soc
weight = ones(n_part,1)/n_part;

%Estimadores
x1_est_mmse = zeros(L,1);
x1_est_mmse(1) = mean(particle(:,1));
x2_est_mmse = zeros(L,1); %soc
x2_est_mmse(1) = mean(particle(:,2));

v_model = zeros(n_part,1);

for k=2:L

    for i=1:n_part
        
        %Modelo de voltaje para la particula anterior (k-1)
        v_model(i) = vl + (v0-vl)*exp(gamma*(particle(i,2)-1))...
            + alfa*vl*(particle(i,2)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(particle(i,2))))...
            - i_meas(k-1)*particle(i,1);
        
        %Importance sampling (prediccion desde k-1 hacia k)
        r1 = sqrt(w_1)*randn;
        r2 = sqrt(w_2)*randn;
        particle_pred(i,1) = particle(i,1) + r1;
        particle_pred(i,2) = particle(i,2) - v_model(i)*i_meas(k-1)*dt/Ecrit + r2;
        if particle_pred(i,2)<0
            particle_pred(i,2) = 0;
        end
        
        %Weight update (medicion de valor en k)
        v_model(i) = vl + (v0-vl)*exp(gamma*(particle_pred(i,2)-1))...
            + alfa*vl*(particle_pred(i,2)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(particle_pred(i,2))))...
            - i_meas(k)*particle_pred(i,1); %prediccion de medicion en k
        innov = v_meas(k) - v_model(i); %innovacion
        %weight(i) = exp( -(( innov )^2)/(2*Rnn) )/(sqrt(2*pi*Rnn));
        weight(i) = exp( -log(sqrt(2*pi*Rnn)) -(( innov )^2)/(2*Rnn) );
    end
    if sum(weight)==0
        %disp(weight);
        disp(innov);
        disp('Error');
        disp(k);
        disp(i);
    end
    %Normalizar pesos
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
            particle(i,:) = particle_pred(ind,:);
        end
    else
        for i=1:n_part
            particle(i,:) = particle_pred(i,:);
        end
        
    end
    
    
    %Estimacion del estado
    x1_est_mmse(k) = mean(particle(:,1));
    x2_est_mmse(k) = mean(particle(:,2));
     
end

figure, plot(t,x2_est_mmse)
hold on
plot(t, x2)
grid, xlabel('Time [s]'), ylabel('SOC'), title('Particle-Filtered SOC ');
legend('Estimate','Ground Truth');


figure
subplot(2,1,1)
plot(t,x2);
hold on
plot(t,x2_est_mmse);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('SOC Estimation with $$x^{(2)}_0=$$ %1.1f and $$P_0^{(2)}=$$ %1.0d',x2_est_mmse(1), 1e-5));
legend('True State','BPF Estimate');


subplot(2,1,2)
error = x2-x2_est_mmse;
plot(t,error);
ylim([-4e-2,4e-2]);
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('SOC Estimation Error');




%%
%Filtrado
%--------------EKF---------------

%Ruidos filtro
w_1 = sigma_1;
w_2 = sigma_2;
Rnn = sigma_n;
Rww = [w_1,0;0,w_2];

%Inicializacion de arreglos
x1_est   = zeros(L,1); %state estimate
x1_pred  = zeros(L,1); %state estimate prediction
P2_est =  zeros(L,1);
x2_est   = zeros(L,1); %state estimate soc
x2_pred  = zeros(L,1); %state estimate prediction soc
v_pred  = zeros(L,1); %measurement prediction
innov   = zeros(L,1); %innovation
R_innov = zeros(L,1); %innovation covariance

%Filter initial conditions
x1_est(1) = Zp;
x2_est(1) = 0.8;
P_est = [0.0001,0;0,0.1];
P2_est(1) = P_est(2,2);

v_model = zeros(n_part,1);


for k=2:L
    
    %Prediction k-1 -> k
    v_model = vl + (v0-vl)*exp(gamma*(x2_est(k-1)-1))...
            + alfa*vl*(x2_est(k-1)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2_est(k-1))))...
            - i_meas(k-1)*x1_est(k-1); %voltaje en k-1
    x1_pred(k) = x1_est(k-1);
    x2_pred(k) = x2_est(k-1) - v_model*i_meas(k-1)*dt/Ecrit;
    [A,C] = jacobians(x1_pred(k),x2_pred(k),i_meas(k-1),gamma, alfa, v0, vl, beta, dt, Ecrit);
    P_pred = A*P_est*(A') + Rww;
    v_pred(k) = vl + (v0-vl)*exp(gamma*(x2_pred(k)-1))...
            + alfa*vl*(x2_pred(k)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2_pred(k))))...
            - i_meas(k)*x1_pred(k); %prediccion de voltaje en k
       
    %Innovation
    innov(k) = v_meas(k) - v_pred(k);
    R_innov(k) = C*P_pred*(C') + Rnn;
    
    %Kalman Gain
    Kgain = P_pred*(C')/R_innov(k);
    
    %Update
    x1_est(k) = x1_pred(k) + Kgain(1)*innov(k);
    x2_est(k) = x2_pred(k) + Kgain(2)*innov(k);
    P_est = ( [1,0;0,1] - Kgain*C )*P_pred;
    P2_est(k) = P_est(2,2);
    

end

figure
subplot(2,1,1)
plot(t,x2);
hold on
plot(t,x2_est);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('SOC Estimation with $$x^{(2)}_0=$$ %1.1f and $$P_0^{(2)}=$$ %1.0d',x2_est(1), P2_est(1)));
legend('True State','EKF Estimate');


subplot(2,1,2)
error = x2-x2_est;
plot(t,error);
hold on
plot(t,2*sqrt(P2_est));
plot(t,-2*sqrt(P2_est));
ylim([-4e-2,4e-2]);
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('SOC Estimation Error');
legend('Error', '$$+2\sqrt{P^{(2)}_{est}} $$','$$-2\sqrt{P^{(2)}_{est}} $$');


figure, plot(t,x2_est)
hold on
plot(t, x2)
grid, xlabel('Time [s]'), ylabel('SOC'), title('Extended Kalman-Filtered SOC ');
legend('Estimate','Ground Truth');
