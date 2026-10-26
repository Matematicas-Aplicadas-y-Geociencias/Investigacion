%% Graficar datos estructurados
clear 
clc
close all

load('Datos_Presion.mat')

%% ── 1. Preparar los datos ────────────────────────────────────────────────────

% Unir todas las fuentes de SF01
time_SF01 = [data_EM.SF01.time.local_1;  ...
             data_EM.SF01.time.local_2;  ...
             data_EM.SF01.time.BigQuery];

PAtm_SF01 = [data_EM.SF01.PAtm.local_1;  ...
             data_EM.SF01.PAtm.local_2;  ...
             data_EM.SF01.PAtm.BigQuery];

% Redondear a 10 minutos — forma simple
time_SF01_r    = dateshift(time_SF01,          'start', 'minute');
time_Conagua_r = dateshift(data_Conagua.time,  'start', 'minute');
time_EC09_r    = dateshift(data_Diverlink.EC09.time, 'start', 'minute');

% Quitar los segundos colgados redondeando al intervalo de 10 min más cercano
time_SF01_r    = dateshift(time_SF01_r,    'start', 'minute', 'nearest');
time_Conagua_r = dateshift(time_Conagua_r, 'start', 'minute', 'nearest');
time_EC09_r    = dateshift(time_EC09_r, 'start', 'minute', 'nearest');

% Ordenar y quitar duplicados
[time_SF01_r, idx] = unique(time_SF01_r);
PAtm_SF01 = PAtm_SF01(idx);

% Grilla regular cada 10 min
t_grid = (min(time_SF01_r) : minutes(10) : max(time_SF01_r))';

% SF01 en la grilla (huecos = NaN)
P_SF01 = NaN(size(t_grid));
[~, ia, ib] = intersect(t_grid, time_SF01_r);
P_SF01(ia) = PAtm_SF01(ib);

% Conagua en la misma grilla
P_Conagua = NaN(size(t_grid));
[~, ia, ib] = intersect(t_grid, time_Conagua_r);
P_Conagua(ia) = data_Conagua.PAtm(ib);

% EC09 en la misma grilla
P_EC09 = NaN(size(t_grid));
[~, ia, ib] = intersect(t_grid, time_EC09_r);
P_EC09(ia) = data_Diverlink.EC09.PAtm(ib);

% Verificar intersección de datos
fprintf('Puntos SF01 con datos:    %d\n', sum(~isnan(P_SF01)))
fprintf('Puntos Conagua con datos: %d\n', sum(~isnan(P_Conagua)))
fprintf('Puntos EC09 con datos: %d\n', sum(~isnan(P_EC09)))
fprintf('Puntos en común:          %d\n', sum(~isnan(P_SF01) & ~isnan(P_Conagua) & ~isnan(P_EC09)))

%% ── 2. Construir dataset de entrenamiento ────────────────────────────────────
% Solo periodos donde ambas tienen datos
mask_train = ~isnan(P_SF01) & (~isnan(P_Conagua) | ~isnan(P_EC09));

P_ref = NaN(size(t_grid));
P_ref(~isnan(P_Conagua)) = P_Conagua(~isnan(P_Conagua));
P_ref(isnan(P_Conagua) & ~isnan(P_EC09)) = P_EC09(isnan(P_Conagua) & ~isnan(P_EC09));

X_train = P_ref(mask_train);   % entrada: Conagua
Y_train = P_SF01(mask_train);      % salida:  SF01

% Normalizar (importante para LSTM)
mu_x = mean(X_train);  sigma_x = std(X_train);
mu_y = mean(Y_train);  sigma_y = std(Y_train);

X_norm = (X_train - mu_x) / sigma_x;
Y_norm = (Y_train - mu_y) / sigma_y;

% Separar 80% entrenamiento / 20% validación
n       = numel(X_norm);
n_train = round(0.8 * n);

% Formato que espera LSTM: cell array de secuencias
XTrain = X_norm(1:n_train)';   % 1 x N
YTrain = Y_norm(1:n_train)';   % 1 x N

XVal   = X_norm(n_train+1:end)';
YVal   = Y_norm(n_train+1:end)';

%% ── 3. Definir la red LSTM ───────────────────────────────────────────────────
layers = [ ...
    sequenceInputLayer(1)
    lstmLayer(64)           % 64 unidades de memoria
    dropoutLayer(0.2)
    fullyConnectedLayer(1)  % salida: 1 valor (presión)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',          150,   ...
    'MiniBatchSize',      512,  ...
    'InitialLearnRate',   1e-3, ...
    'GradientThreshold',  1,    ...
    'ValidationData',     {XVal, YVal}, ...
    'ValidationFrequency', 10, ...
    'ExecutionEnvironment','gpu', ...  % ← aquí usa la GPU
    'Verbose',            false, ...
    'Plots',              'training-progress');

%% ── 4. Entrenar ─────────────────────────────────────────────────────────────
net = trainNetwork(XTrain, YTrain, layers, options);

%% ── 5. Rellenar los huecos ───────────────────────────────────────────────────
mask_hueco = isnan(P_SF01) & ~isnan(P_ref);
fprintf('Huecos a rellenar: %d\n', sum(mask_hueco))

X_hueco = (P_ref(mask_hueco) - mu_x) / sigma_x;
Y_pred_norm = predict(net, {X_hueco'});
Y_pred = Y_pred_norm{1}' * sigma_y + mu_y;   % desnormalizar

P_SF01_relleno = P_SF01;
P_SF01_relleno(mask_hueco) = Y_pred;

%% ── 6. Graficar ──────────────────────────────────────────────────────────────
figure
plot(t_grid, P_SF01,                        'k-', 'MarkerSize', 2, 'DisplayName', 'SF01 original')
hold on
plot(t_grid(mask_hueco), P_SF01_relleno(mask_hueco),            'r-', 'MarkerSize', 2, 'DisplayName', 'LSTM relleno')
plot(t_grid, P_Conagua,                     'b-', 'MarkerSize', 2, 'DisplayName', 'Conagua')
plot(t_grid, P_EC09,                                     'g-', 'MarkerSize', 2, 'DisplayName', 'EC09')
legend('show', 'Location', 'best')
title('SF01 — Relleno con LSTM')
ylabel('Presión (hPa)')