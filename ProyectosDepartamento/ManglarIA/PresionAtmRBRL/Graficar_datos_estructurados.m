%% Graficar datos estructurados
clear 
close all

load('Datos_Presion.mat')

%% primero obtengo el número de estaciones

estacion_names = fieldnames(data_EM);
num_estaciones = numel(estacion_names);

for i=1:num_estaciones


%% Primero grafico los datos de la CONAGUA

figure(i)
plot(data_Conagua.time, data_Conagua.PAtm, 'k.', 'DisplayName', 'EMA Conagua')
hold on

%% y los datos de EC09 (u otra referencia conveniente).
plot(data_Diverlink.EC09.time, data_Diverlink.EC09.PAtm, 'r.', 'DisplayName', 'DivLink EC09')

%% Ahora voy a obtener el número de datos de BigQuery (debería ser solo 1, pero lo dejo para ir actualizando solo con lo "extra" en el futuro)

    p_fields = fieldnames(data_EM.(estacion_names{i}).PAtm);
    starts_BigQuery = strncmpi(p_fields,'BigQu',5); % compara primeras 5 letras, case-insensitive
    num_BigQ = sum(starts_BigQuery);

if num_BigQ>0
plot(data_EM.(estacion_names{i}).time.BigQuery, data_EM.(estacion_names{i}).PAtm.BigQuery, 'b.', 'DisplayName','BigQuery')
end 

%% ahora voy a obtener el número de datos locales que hay en la estación i
     p_fields = fieldnames(data_EM.(estacion_names{i}).PAtm);
    starts_local = strncmpi(p_fields,'local',5); % compara primeras 5 letras, case-insensitive
    num_local = sum(starts_local);

%% y también el número de archivos de Weatherlink (Wlink), que por ahora son 2, pero lo dejo asi por si se consiguen más 
    starts_Wlink = strncmpi(p_fields,'Wlink',5); % compara primeras 5 letras, case-insensitive
    num_Wlink = sum(starts_Wlink);

%% ahora toca graficar los archivos locales

for k=1: num_local
% plot(data_EM.(ID_esta))
name=['local_', num2str(k)];
plot(data_EM.(estacion_names{i}).time.(name), data_EM.(estacion_names{i}).PAtm.(name), 'g.', 'DisplayName','Local')

end 
   
%%  y ahora los archivos con WLink

if num_Wlink>=1
    
    for j=1:num_Wlink
    name=['Wlink_', num2str(j)];
    plot(data_EM.(estacion_names{i}).time.(name), data_EM.(estacion_names{i}).PAtm.(name), 'm.', 'DisplayName','WeatherLink')

    end 

end
%Primero 



title(estacion_names{i})
legend('show','Location','best')
set(gca,'FontSize',16)

% Eliminar entradas de leyenda repetidas
legh = findobj(gcf, '-property', 'DisplayName');           % objetos con DisplayName
names = get(legh, 'DisplayName');                         % obtener nombres
if iscell(names)
    [~, ia] = unique(names, 'stable');                    % índices únicos (mantener orden)
    keep = legh(ia);
else
    keep = legh;                                          % un único objeto
end
% Crear leyenda solo con objetos únicos
if ~isempty(keep)
    legend(keep, get(keep, 'DisplayName'), 'Location', 'best')
end

% Guardar la figura en .fig y .jpg con el nombre requerido
fname = ['RBRL-', estacion_names{i}, '_PAtm'];
savefig(gcf, [fname, '.fig'])
print(gcf, [fname, '.jpg'], '-djpeg', '-r300')


end

