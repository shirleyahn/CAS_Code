%% Penta Alanine Spectral Clustering
load('ball_clustering_15.txt')
ball_clustering = ball_clustering_15;
m = size(ball_clustering, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
  a(i) = ball_clustering(i, 5);
  plot(ball_clustering(i, 1), ball_clustering(i, 2),'ro', 'MarkerSize', 1.0);
end
s = 0.2; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(ball_clustering(:,1), ball_clustering(:,2), 1, a, 'Linewidth', 1.0), colormap(jet), colorbar;
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\phi (deg)');
ylabel('\psi (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Penta Alanine (first)
load('total_weight_on_each_ball_65.txt')
total_weight = total_weight_on_each_ball_65;
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 1), total_weight(i, 2), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 1), total_weight(:, 2), 1, -0.0019872041*300*log(a), 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 25]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\phi (deg)');
ylabel('\psi (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Penta Alanine (second)
load('total_weight_on_each_ball_65.txt')
total_weight = total_weight_on_each_ball_65;
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 3), total_weight(i, 4), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 3), total_weight(:, 4), 1, -0.0019872041*300*log(a), 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 25]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\phi (deg)');
ylabel('\psi (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Penta Alanine (third)
load('total_weight_on_each_ball_65.txt')
total_weight = total_weight_on_each_ball_65;
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 5), total_weight(i, 6), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 5), total_weight(:, 6), 1, -0.0019872041*300*log(a), 'Linewidth', 1.0), colormap(jet), colorbar, caxis([0 25]);
axis([-180.0 180.0 -180.0 180.0]);
xlabel('\phi (deg)');
ylabel('\psi (deg)');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);
