%% Triazine Diffusion Map
diffusion_map_coords = load('diffusion_map_coords_64.txt');
m = size(diffusion_map_coords, 1);
a = zeros(m, 1);
for i=1:m
  a(i) = diffusion_map_coords(i, 5);
end
figure;scatter(diffusion_map_coords(:, 2), diffusion_map_coords(:, 3), 1000, -0.0019872041*300*log(a)), colormap(jet), colorbar;
figure;scatter3(diffusion_map_coords(:, 2), diffusion_map_coords(:, 3), diffusion_map_coords(:, 4), 1000, -0.0019872041*300*log(a)), colormap(jet), colorbar;


%% Triazine (first)
total_weight = load('total_weights_on_each_bin_142.txt');
m=size(total_weight, 1);
a=zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 5);
    plot(total_weight(i, 1), total_weight(i, 2), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 1), total_weight(:, 2), 1, -0.0019872041*300*log(a), 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 100]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\omega_1 (deg)');
ylabel('\omega_2 (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);

%% Triazine (first)
total_weight = load('total_weights_on_each_bin_142.txt');
m=size(total_weight, 1);
a=zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 5);
    plot(total_weight(i, 3), total_weight(i, 4), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 3), total_weight(:, 4), 1, -0.0019872041*300*log(a), 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 100]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\omega_1 (deg)');
ylabel('\omega_2 (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Triazine Initial Conditions
x = (-2.332);
y = (-4.755);
a = (0.0);
figure; hold on
plot(x, y, 'ro', 'MarkerSize', 1.0);
s = 20.0; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(x, y, 1, a, 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 100]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\omega_1 (deg)');
ylabel('\omega_2 (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);
