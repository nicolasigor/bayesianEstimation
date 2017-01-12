%Leer datos del experimento del paper, train y vali de la bateria 1

set(0,'defaulttextinterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

train1 = csvread('2015_train1_current.csv');
train2 = csvread('2015_train1_voltage.csv');
test1 = csvread('2015_vali1_current.csv');
test2 = csvread('2015_vali1_voltage.csv');

%Extraer series de datos
i_train = train1(:,2);
v_train = train2(:,2);
L_train = length(i_train);
i_test = test1(:,2);
v_test = test2(:,2);
v_test = v_test-0.1; %correction due to the extraction of the original data
L_test = length(i_test);

%Crear arreglos de tiempo que comienzan de 0
dt = 2; %paso de tiempo en segundos
t_train = dt*(0:(L_train-1));
t_test = dt*(0:(L_test-1));

%Mostrar los datos
figure %Entrenamiento
subplot(1,2,1)
plot(t_train,i_train);
grid, xlabel('Time [s]'), ylabel('Current [A]'), title('Battery \#1 Measured Discharge Current');
subplot(1,2,2)
plot(t_train,v_train);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \#1 Measured Discharge Voltage');

figure %Validacion
subplot(1,2,1)
plot(t_test,i_test);
grid, xlabel('Time [s]'), ylabel('Current [A]'), title('Battery \#1 Measured Discharge Current');
subplot(1,2,2)
plot(t_test,v_test);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \#1 Measured Discharge Voltage');

%%
%Estimar parametros.
%Parametros estimados en el paper para la bateria 1
%alfa = 0.08, beta =16, gamma = 19.65,  v0 = 4.12, vl = 3.987, Ecrit = 20127, Zp = 0.30

%------Estimacion de Ecrit--------

%Mostrar energia entregada por la bateria en el entrenamiento
energyDel = [];
for k=1:L_train
    energyDel = [energyDel, sum(i_train(1:k).*v_train(1:k)*dt) ];
end
% figure
% plot(t_train,energyDel);
% grid, xlabel('Time [s]'), ylabel('Energy [J]'), title('Battery #1 Measured Energy Delivered');

Ecrit = energyDel(end); %20119
fprintf('Ecrit: %f \n',Ecrit); 

%%
%------Estimacion de Zp--------
%en los tiempos 740, 800 el segundo del primer par de pulsos
di1 = 2.399 - 1.401;
dv1 = 3.196 - 3.464;
zp1 = abs(dv1/di1);
%en los tiempos 1770, 1830 el segundo del segundo par de pulsos
di2 = 2.395 - 1.404;
dv2 = 3.077 - 3.356;
zp2 = abs(dv2/di2);

Zp = mean([zp1,zp2]); %0.275
fprintf('Zp: %f \n',Zp); 

%%
%------Arreglos OCV vs SOC --------
soc = [1];
for k = 2:L_train
    soc = [soc; soc(k-1) - v_train(k)*i_train(k)*dt/Ecrit];
end

ocv = v_train + i_train*Zp;
% figure
% plot(soc,ocv);
% grid, xlabel('SOC'), ylabel('Voltage [V]'), title('Battery #1 Discharge Open Circuit Voltage');
% set(gca,'xdir','reverse')

%------Encontrar zonas --------
zone1soc = soc( soc>0.7 );
zone2soc = soc( soc<=0.7 & soc>0.25 );
zone3soc = soc( soc<=0.25 );
%ind1end = length(zone1);
%ind2end = L_train - length(zone3);
zone1ocv = ocv( soc>0.7 );
zone2ocv = ocv( soc<=0.7 & soc>0.25 );
zone3ocv = ocv( soc<=0.25 );

figure, plot(zone1soc,zone1ocv,'LineWidth',1);
hold on
plot(zone2soc,zone2ocv,'LineWidth',1);
plot(zone3soc,zone3ocv,'LineWidth',1);
grid, xlabel('SOC'), ylabel('Voltage [V]'), title('Battery \#1 Discharge Open Circuit Voltage');
set(gca,'xdir','reverse')
legend('Zone 1','Zone 2','Zone 3')

%%
%------Estimacion de vl y alfa (zona 2) --------

coef = polyfit(zone2soc, zone2ocv, 1);
a = coef(1);
b = coef(2);
alfa = 1/(1+(b/a)); % 0.067
vl = a/alfa; % 3.987
fprintf('alfa: %f \n',alfa); 
fprintf('v_L: %f \n',vl); 

%%
%------Estimacion de v0 y gamma (zona 1) --------

v0 = ocv(1); % 4.0435
notGamma = vl + alfa*vl*(zone1soc-1) - zone1ocv;
zone1rms = @(x) sqrt(sum(( notGamma + (v0-vl)*exp(x*(zone1soc-1)) ).^2));

gammaAux = linspace(0,30,200);
rmsGammaAux = [];
for k = 1:200
    rmsGammaAux = [rmsGammaAux, zone1rms(gammaAux(k))];
end
figure, plot(gammaAux, rmsGammaAux);
grid, xlabel('gamma'), ylabel('RMS'), title('Model RMS vs gamma in zone 1');

gamma = fminsearch(zone1rms,15); %encontrar gamma 14.825
fprintf('v_0: %f \n',v0);
fprintf('gamma: %f \n',gamma); 


%%
%------Estimacion de beta (zona 3) --------

notBeta = vl + (v0-vl)*exp(gamma*(zone3soc-1)) + alfa*vl*(zone3soc-1) - zone3ocv;
zone3rms = @(x) sqrt(sum(( notBeta + (1-alfa)*vl*( exp(-x)-exp(-x*sqrt(zone3soc)) ) ).^2));

betaAux = linspace(0,30,200);
rmsBetaAux = [];
for k = 1:200
    rmsBetaAux = [rmsBetaAux, zone3rms(betaAux(k))];
end
figure, plot(betaAux, rmsBetaAux);
grid, xlabel('beta'), ylabel('RMS'), title('Model RMS vs beta in zone 3');

beta = fminsearch(zone3rms,15); %encontrar beta 15.026
fprintf('beta: %f \n',beta);

%%
% gamma and beta
figure
subplot(1,2,1)
plot(gammaAux, rmsGammaAux);
grid, xlabel('$$\gamma$$'), ylabel('RSE'), title('Model RSE vs $$\gamma$$ in zone 1');
subplot(1,2,2)
plot(betaAux, rmsBetaAux);
grid, xlabel('$$\beta$$'), ylabel('RSE'), title('Model RSE vs $$\beta$$ in zone 3');

%%
%Resumen de parametros
fprintf('\n');
disp('Parameters of the model vs Pola parameters:');
fprintf('alfa: %f (0.08)\n',alfa);
fprintf('beta: %f (16)\n',beta);
fprintf('gamma: %f (19.65)\n',gamma); 
fprintf('v_0: %f (4.12)\n',v0);
fprintf('v_L: %f (3.987)\n',vl); 
fprintf('Ecrit: %f (20127)\n',Ecrit); 
fprintf('Zp: %f (0.30)\n',Zp); 


%%
%Contrastar modelo con los datos del entrenamiento OCV vs SOC
ocv_model = vl + (v0-vl)*exp(gamma*(soc-1)) + alfa*vl*(soc-1)...
    + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(soc)));
figure
plot(soc,ocv);
grid, xlabel('SOC'), ylabel('Voltage [V]'), title('Battery #1 Discharge Open Circuit Voltage');
set(gca,'xdir','reverse')
hold on
plot(soc,ocv_model);
legend('Measured','Model');

%%
%Contrastar modelo con los datos de voltaje medidos sin ruido

v_train_model = zeros(L_train,1);
x1 = zeros(L_train,1);
x2 = zeros(L_train,1);
%Initial conditions
x1(1) = Zp;
x2(1) = 1;
v_train_model(1) = v0 - i_train(1)*x1(1);
for k = 2:L_train
    %State transition
    x1(k) = x1(k-1);
    x2(k) = x2(k-1) - v_train_model(k-1)*i_train(k-1)*dt/Ecrit;
    %Measurement
    v_train_model(k) = vl + (v0-vl)*exp(gamma*(x2(k)-1))...
        + alfa*vl*(x2(k)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2(k))))...
        - i_train(k)*x1(k);
end

figure
subplot(2,1,1)
plot(t_train,v_train);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \# 1 Measured Discharge Voltage');
hold on
plot(t_train,v_train_model);
legend('Measured','Model');

subplot(2,1,2)
plot(t_train,v_train-v_train_model);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \# 1 Error between model and measurements in voltage');

%mesStd = std(v_train-v_train_model);

%%
%Contrastar modelo con los datos de voltaje medidos sin ruido de test

v_test_model = zeros(L_test,1);
x1 = zeros(L_test,1);
x2 = zeros(L_test,1);
%Initial conditions
x1(1) = Zp;
x2(1) = 1;
v_test_model(1) = v0 - i_test(1)*x1(1);
%vl = 3.987;
%Ecrit = 20119+700;
for k = 2:L_test
    %State transition
    x1(k) = x1(k-1);
    x2(k) = x2(k-1) - v_test_model(k-1)*i_test(k-1)*dt/(Ecrit);
    %Measurement
    v_test_model(k) = vl + (v0-vl)*exp(gamma*(x2(k)-1))...
        + alfa*vl*(x2(k)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2(k))))...
        - i_test(k)*x1(k);
end

figure
subplot(2,1,1)
plot(t_test,v_test);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \# 1 Measured Discharge Voltage Test');
hold on
plot(t_test,v_test_model);
legend('Measured','Model');


subplot(2,1,2)
plot(t_test,v_test-v_test_model);
grid, xlabel('Time [s]'), ylabel('Voltage [V]'), title('Battery \# 1 Error between model and measurements in voltage');

%%
%Validation

%GroundTruth offline
EcritTest = sum(i_test.*v_test*dt);
socTest = zeros(L_test,1);
socTest(1) = 1;
for k = 2:L_test
    socTest(k) = socTest(k-1) - v_test(k)*i_test(k)*dt/EcritTest;
end

%%
%--------------PF Bootstrap---------------
n_part = 100;

%Ruidos
sigma_1 = 1e-8;
sigma_2 = 1e-5;
Rnn = 1e-2;

%Inicializacion de particulas
particle = zeros(n_part,2);
particle_pred =  zeros(n_part,2);
particle(:,1) = ones(n_part,1)*Zp;
particle(:,2) = unifrnd(0.8,0.9, n_part, 1);
weight = ones(n_part,1)/n_part;

%Estimadores
%x1_est_map = zeros(L_test,1);
x1_est_mmse = zeros(L_test,1);
x1_est_mmse(1) = mean(particle(:,1));
%x2_est_map = zeros(L_test,1); %soc
x2_est_mmse = zeros(L_test,1); %soc
x2_est_mmse(1) = mean(particle(:,2));

v_model = zeros(n_part,1);

for k=2:L_test
    if k == 350
        sigma_2 = 1e-6;
    end
    if k == 450
        sigma_2 = 1e-7;
    end
    if k == 550
        sigma_2 = 1e-8;
    end
    
    for i=1:n_part
        
        %Modelo de voltaje para la particula anterior (k-1)
        v_model(i) = vl + (v0-vl)*exp(gamma*(particle(i,2)-1))...
            + alfa*vl*(particle(i,2)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(particle(i,2))))...
            - i_test(k-1)*particle(i,1);
        
        %Importance sampling (prediccion desde k-1 hacia k)
        w1 = sqrt(sigma_1)*randn;
        w2 = sqrt(sigma_2)*randn;
        particle_pred(i,1) = particle(i,1) + w1;
        particle_pred(i,2) = particle(i,2) - v_model(i)*i_test(k-1)*dt/Ecrit + w2;
        if particle_pred(i,2)<0
            particle_pred(i,2) = 0;
        end
        
        %Weight update (medicion de valor en k)
        v_model(i) = vl + (v0-vl)*exp(gamma*(particle_pred(i,2)-1))...
            + alfa*vl*(particle_pred(i,2)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(particle_pred(i,2))))...
            - i_test(k)*particle_pred(i,1); %prediccion de medicion en k
        innov = v_test(k) - v_model(i); %innovacion
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
    if N_eff < n_part
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
    %temp1 = particle( (weight==max(weight)) ,1);
    %temp2 = particle( (weight==max(weight)) ,2);
    %x1_est_map(k) = temp1(1);
    %x2_est_map(k) = temp2(1);
     
end

figure, plot(t_test,x2_est_mmse)
hold on
plot(t_test, x2)
grid, xlabel('Time [s]'), ylabel('SOC'), title('Particle-Filtered SOC ');
legend('Estimate','Ground Truth');


figure
subplot(2,1,1)
plot(t_test,x2);
hold on
plot(t_test,x2_est_mmse);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('SOC Estimation with $$x^{(2)}_0=$$ %1.1f and $$P_0^{(2)}=$$ %1.0d,%1.0d',x2_est_mmse(1), 1e-5,1e-8));
legend('True State','BPF Estimate');


subplot(2,1,2)
error = x2-x2_est_mmse;
plot(t_test,error);
ylim([-4e-2,4e-2]);
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('SOC Estimation Error');

%%

%Filtrado
%--------------EKF---------------

%Noise
sigma_1 = 1e-8;
sigma_2 = 1e-8;
sigma_n = 1e-4;

%Ruidos filtro
w_1 = sigma_1;
w_2 = sigma_2;
Rnn = sigma_n;
Rww = [w_1,0;0,w_2];

%Inicializacion de arreglos
x1_est   = zeros(L_test,1); %state estimate
x1_pred  = zeros(L_test,1); %state estimate prediction
P2_est =  zeros(L_test,1);
x2_est   = zeros(L_test,1); %state estimate soc
x2_pred  = zeros(L_test,1); %state estimate prediction soc
v_pred  = zeros(L_test,1); %measurement prediction
innov   = zeros(L_test,1); %innovation
R_innov = zeros(L_test,1); %innovation covariance

%Filter initial conditions
x1_est(1) = Zp;
x2_est(1) = 0.8;
P_est = [0.0001,0;0,0.1];
P2_est(1) = P_est(2,2);

v_model = zeros(n_part,1);


for k=2:L_test
    
    %Prediction k-1 -> k
    v_model = vl + (v0-vl)*exp(gamma*(x2_est(k-1)-1))...
            + alfa*vl*(x2_est(k-1)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2_est(k-1))))...
            - i_test(k-1)*x1_est(k-1); %voltaje en k-1
    x1_pred(k) = x1_est(k-1);
    x2_pred(k) = x2_est(k-1) - v_model*i_test(k-1)*dt/Ecrit;
    if x2_pred(k)<0
            x2_pred(k) = 0;
    end
    [A,C] = jacobians(x1_pred(k),x2_pred(k),i_test(k-1),gamma, alfa, v0, vl, beta, dt, Ecrit);
    P_pred = A*P_est*(A') + Rww;
    v_pred(k) = vl + (v0-vl)*exp(gamma*(x2_pred(k)-1))...
            + alfa*vl*(x2_pred(k)-1) + (1-alfa)*vl*( exp(-beta)-exp(-beta*sqrt(x2_pred(k))))...
            - i_test(k)*x1_pred(k); %prediccion de voltaje en k
       
    %Innovation
    innov(k) = v_test(k) - v_pred(k);
    R_innov(k) = C*P_pred*(C') + Rnn;
    
    %Kalman Gain
    Kgain = P_pred*(C')/R_innov(k);
    
    %Update
    x1_est(k) = x1_pred(k) + Kgain(1)*innov(k);
    x2_est(k) = x2_pred(k) + Kgain(2)*innov(k);
    if x2_est(k)<0
        x2_est(k) = 0;
    end
    P_est = ( [1,0;0,1] - Kgain*C )*P_pred;
    P2_est(k) = P_est(2,2);
    

end

figure
subplot(2,1,1)
plot(t,socTest);
hold on
plot(t,x2_est);
grid, xlabel('Time [s]'), ylabel('State'), title(sprintf('SOC Estimation with $$x^{(2)}_0=$$ %1.1f and $$P_0^{(2)}=$$ %1.0d',x2_est(1), P2_est(1)));
legend('True State','EKF Estimate');


subplot(2,1,2)
error = socTest-x2_est;
plot(t,error);
hold on
plot(t,2*sqrt(P2_est));
plot(t,-2*sqrt(P2_est));
ylim([-5e-2,5e-2]);
grid, xlabel('Time [s]'), ylabel('Error');
%title(sprintf('State Estimation Error with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x_est(1), P_est(1)));
title('SOC Estimation Error');
legend('Error', '$$+2\sqrt{P^{(2)}_{est}} $$','$$-2\sqrt{P^{(2)}_{est}} $$');


figure, plot(t,x2_est)
hold on
plot(t, socTest)
grid, xlabel('Time [s]'), ylabel('SOC'), title('Extended Kalman-Filtered SOC ');
legend('Estimate','Ground Truth');




