%% MakeDemoPlots
% A script to plot simple Ising model simulation results.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License

%% Evaluate the theoretical energies and magnetizations at each temperature
% (see http://www.scholarpedia.org/article/Ising_model:_exact_results)

% Compute the exact value of the critical temperature for the infinite
% lattice
Tc = 2 / log(1 + sqrt(2));

% Define an anonymous function to compute the exact energy for the infinite
% lattice
kappa = @(T) 2 * sinh( 2 ./ T ) ./ (cosh( 2 ./ T ).^2);
K = @(T) integral( @(phi) 1./ sqrt( 1 - (kappa(T).^2) .* (sin(phi).^2) ), 0, pi/2 , 'ArrayValued', true);
E_theoretical = @(T) -2 * tanh(2./T) - ((sinh(2./T).^2 - 1) ./ (sinh(2./T).*cosh(2./T))) .* ( (2/pi) * K(T) - 1);

% Define an anonymous function to compute the exact magnetization for the infinite
% lattice
M_theoretical = @(T) ((1 - (sinh(2./T)).^(-4)).^(1/8)) .* ((T - Tc) < 0);

%% Compute statistics of simulations

% Mean energy
E = mean(E_iter(nBurnin:end,:));

% Standard deviation of energy
dE = std(abs(E_iter(nBurnin:end,:))) / sqrt(size(E,1));

% Difference between empirical + theoretical energies
E_residual = (E - E_theoretical(T));

% Mean magnitization
M = mean(abs(M_iter(nBurnin:end,:)));

% Standard deviation of magnetization
dM = std(abs(M_iter(nBurnin:end,:))) / sqrt(size(E,1));

% Difference between empirical + theoretical magnetizations
M_residual = (M - M_theoretical(T));

% Specific heat
cv = var(E_iter(nBurnin:end,:))./(T.^2);

% Magnetic susceptibility
chi = var(abs(M_iter(nBurnin:end,:))) ./ T;

%% Plot everything

% Define query points at which to plot theoretical values
x = linspace(min(T), max(T), length(T)*10);

% Set color order
colorOrder = lines(3);

% Plot the measured and theoretical energies as a function of temperature
f = figure('Position',[200,500,1000,1000],'WindowStyle','docked');
plot(T, E, 'linewidth', 2, 'Color', colorOrder(1,:));
hold on;
p = patch([T, fliplr(T)], [E + dE, fliplr(E-dE)], 1, 'FaceColor',colorOrder(1,:),'LineStyle','None', 'FaceAlpha', 0.25);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(x, E_theoretical(x), 'linewidth', 2, 'Color', colorOrder(2,:));
axis('square','tight');
xlabel('T (J/k_B)');
ylabel('E (J)');
set(gca,'FontSize',20);
set(gca,'LineWidth', 2);
set(gca,'Box', 'off');
legend({'measured','theoretical'}, 'location','northwest');

% Plot the residuals in an inset
inset = axes(f, 'Position', [0.4 0.2 0.4 0.4]);
plot(T, E_residual, 'linewidth', 2, 'Color', colorOrder(3,:));
hold on;
patch([T, fliplr(T)], [E_residual + dE, fliplr(E_residual-dE)], 1, 'FaceColor',colorOrder(3,:),'LineStyle','None', 'FaceAlpha', 0.25);
xlabel('T (J/k_B)');
ylabel('residual');
axis('square');
set(inset, 'FontSize', 15);
set(inset, 'LineWidth', 2);
inset.XAxis.Label.FontSize = 20;
inset.YAxis.Label.FontSize = 20;
xlim([min(x) max(x)]);
set(inset, 'box','off');

% Plot the measured and theoretical mangeitzations as a function of temperature
f = figure('Position',[200,500,1000,1000],'WindowStyle','docked');
plot(T, M, 'linewidth', 2, 'Color', colorOrder(1,:));
hold on;
p = patch([T, fliplr(T)], [M + dM, fliplr(M-dM)], 1, 'FaceColor',colorOrder(1,:),'LineStyle','None', 'FaceAlpha', 0.25);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(x, M_theoretical(x), 'linewidth', 2, 'Color', colorOrder(2,:));
axis('square','tight');
xlabel('T (J/k_B)');
ylabel('|M|');
set(gca,'FontSize',20);
set(gca,'LineWidth', 2);
set(gca,'Box', 'off');
legend({'measured','theoretical'}, 'location','southeast');

% Plot the residuals in an inset
inset = axes(f, 'Position', [0.4 0.5 0.4 0.4]);
plot(T, M_residual, 'linewidth', 2, 'Color', colorOrder(3,:));
hold on;
patch([T, fliplr(T)], [M_residual + dM, fliplr(M_residual-dM)], 1, 'FaceColor',colorOrder(3,:),'LineStyle','None', 'FaceAlpha', 0.25);
xlabel('T (J/k_B)');
ylabel('residual');
axis('square');
set(inset, 'FontSize', 15);
set(inset, 'LineWidth', 2);
inset.XAxis.Label.FontSize = 20;
inset.YAxis.Label.FontSize = 20;
xlim([min(x) max(x)]);
set(inset, 'box','off');

% Plot the specific heat as a function of temperature
figure('Position',[200,500,1000,1000],'WindowStyle','docked');
plot(T, cv,  '-x', 'linewidth', 2);
xlabel('T (J/k_B)');
ylabel('c_V');
axis('square');
ConfAxis;

% Plot the magnetic susceptibility as a function of temperature
figure('Position',[200,500,1000,1000],'WindowStyle','docked');
plot(T,  chi,  '-x', 'linewidth', 2);
xlabel('T (J/k_B)');
ylabel('\chi');
axis('square');
ConfAxis;
