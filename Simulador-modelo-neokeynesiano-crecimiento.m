%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%% w= 1 y r = 0.05  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

% Parámetros de la Economía
    w = 1;                       % Elasticidad del capital en la producción
    beta  = 0.99;                       % Factor de descuento del hogar
    r = 0.01;                        % Depreciación del capital

% Parámetros Computacionales
    d     = Inf;                        % Inicializar Fv - V
    tol   = 0.0000000000001;            % Tolerancia
    na    = 1000;                       % Número de puntos en grilla
    a    = linspace(0,10,na);          % Grilla de capital

% % % ========================================== % % %
% % %  iteración
% % % ========================================== % % %

% Se crea matriz para el consumo a partir de valores que toma el capital

    for j = 1:na;
        c(:,j) = w - a(j) +(1+r).*a;
    end
 c = c.*(c >= 0) + eps;

% matriz de utilidad
u = -1 ./ c;
v = zeros(na,1)';                       % Conjetura inicial con ceros

tic;
while (d > tol)

    % Ecuación de Bellman
    RHS = u + beta * v;

    % Se maximiza sobre RHS para obtener Tv y su ubicación
    [Tv, argmax] = max(RHS');

    % Función de política
    g = a(argmax);

    % Distancia entre la iteración actual y la iteración anterior
    d = max(abs(Tv - v));

    % Se actualiza v hasta que converja
    v = Tv;

end
toc

% % % ========================================== % % %
% % %  Gráficos
% % % ========================================== % % %

% Gráfica 1: Función valor.

figure(1);
set(gcf,'color', 'w')
plot(a(2:na), v(2:na),'r','LineWidth', 1)
box on
title('Función valor (w = 1 , r = 0.05)')
legend boxoff
legend('Función Valor','fontsize', 10,'location','best')
xlabel('$k$','Interpreter','latex','fontsize', 12)
ylabel('$\upsilon(k)$','Interpreter','latex','fontsize', 12)

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig1

% Grafica 2: Capital.

figure(2);
set(gcf,'color', 'w')

recta45 = a; 

hold on
plot(a(2:na), g(2:na),'r','LineWidth', 1)
plot(a(2:na), recta45(2:na),'y--','LineWidth', 1)
hold off
box on
title('Capital (w = 1 , r = 0.05)')
legend boxoff
legend('$k_{t+1}$', 'Recta de $45^{\circ}$','fontsize', 10,'Interpreter','latex','location','best')
xlabel('$k_{t}$','fontsize', 20,'Interpreter','latex')
ylabel('$k_{t+1}$','fontsize', 20,'Interpreter','latex')

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%% w= 2 y r = 0.05  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Parámetros de la Economía
    w = 2;                       % Elasticidad del capital en la producción
    beta  = 0.99;                       % Factor de descuento del hogar
    r = 0.05;                        % Depreciación del capital

% Parámetros Computacionales
    d     = Inf;                        % Inicializar Fv - V
    tol   = 0.0000000000001;            % Tolerancia
    na    = 1000;                       % Número de puntos en grilla
    a    = linspace(0,10,na);          % Grilla de capital

% % % ========================================== % % %
% % %  iteración
% % % ========================================== % % %

% Se crea matriz para el consumo a partir de valores que toma el capital

    for j = 1:na
        c(:,j) = w - a(j) +(1+r).*a;
    end
 c = c.*(c >= 0) + eps;

% matriz de utilidad
u = -1 ./ c;
v = zeros(na,1)';                       % Conjetura inicial con ceros

tic;
while (d > tol)

    % Ecuación de Bellman
    RHS = u + beta * v;

    % Se maximiza sobre RHS para obtener Tv y su ubicación
    [Tv, argmax] = max(RHS');

    % Función de política
    g = a(argmax);

    % Distancia entre la iteración actual y la iteración anterior
    d = max(abs(Tv - v));

    % Se actualiza v hasta que converja
    v = Tv;

end
toc

% % % ========================================== % % %
% % %  Gráficos
% % % ========================================== % % %

% Gráfica 1: Función valor.

figure(3);
set(gcf,'color', 'w')
plot(a(2:na), v(2:na),'r','LineWidth', 1)
box on
title('Función valor (w = 2 , r = 0.05)')
legend boxoff
legend('Función Valor','fontsize', 10,'location','best')
xlabel('$k$','Interpreter','latex','fontsize', 12)
ylabel('$\upsilon(k)$','Interpreter','latex','fontsize', 12)

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig3

% Grafica 2: Capital.

figure(4);
set(gcf,'color', 'w')

recta45 = a; 

hold on
plot(a(2:na), g(2:na),'r','LineWidth', 1)
plot(a(2:na), recta45(2:na),'y--','LineWidth', 1)
hold off
box on
title('Capital (w = 2 , r = 0.05)')
legend boxoff
legend('$k_{t+1}$', 'Recta de $45^{\circ}$','fontsize', 10,'Interpreter','latex','location','best')
xlabel('$k_{t}$','fontsize', 20,'Interpreter','latex')
ylabel('$k_{t+1}$','fontsize', 20,'Interpreter','latex')

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%% w= 1 y r = 0.1  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
tic;
% Parámetros de la Economía
    w = 1;                       % Elasticidad del capital en la producción
    beta  = 0.99;                       % Factor de descuento del hogar
    r = 0.1;                        % Depreciación del capital

% Parámetros Computacionales
    d     = Inf;                        % Inicializar Fv - V
    tol   = 0.0000000000001;            % Tolerancia
    na    = 1000;                       % Número de puntos en grilla
    a    = linspace(0,10,na);          % Grilla de capital

% % % ========================================== % % %
% % %  iteración
% % % ========================================== % % %

% Se crea matriz para el consumo a partir de valores que toma el capital

    for j = 1:na
        c(:,j) = w - a(j) +(1+r).*a;
    end
 c = c.*(c >= 0) + eps;

% matriz de utilidad
u = -1 ./ c
v = zeros(na,1)';                       % Conjetura inicial con ceros


while (d > tol)

    % Ecuación de Bellman
    RHS = u + beta * v;

    % Se maximiza sobre RHS para obtener Tv y su ubicación
    [Tv, argmax] = max(RHS');

    % Función de política
    g = a(argmax);

    % Distancia entre la iteración actual y la iteración anterior
    d = max(abs(Tv - v));

    % Se actualiza v hasta que converja
    v = Tv;

end


% % % ========================================== % % %
% % %  Gráficos
% % % ========================================== % % %

% Gráfica 1: Función valor.

figure(5);
set(gcf,'color', 'w')
plot(a(2:na), v(2:na),'r','LineWidth', 1)
box on
title('Función valor (w = 1 , r = 0.1)')
legend boxoff
legend('Función Valor','fontsize', 10,'location','best')
xlabel('$k$','Interpreter','latex','fontsize', 12)
ylabel('$\upsilon(k)$','Interpreter','latex','fontsize', 12)

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig5

% Grafica 2: Capital.

figure(6);
set(gcf,'color', 'w')

recta45 = a; 

hold on
plot(a(2:na), g(2:na),'r','LineWidth', 1)
plot(a(2:na), recta45(2:na),'y--','LineWidth', 1)
hold off
box on
title('Capital (w = 1 , r = 0.1)')
legend boxoff
legend('$k_{t+1}$', 'Recta de $45^{\circ}$','fontsize', 10,'Interpreter','latex','location','best')
xlabel('$k_{t}$','fontsize', 20,'Interpreter','latex')
ylabel('$k_{t+1}$','fontsize', 20,'Interpreter','latex')

orient landscape
print -dpdf -fillpage IteracionFuncionValor_fig6
toc
