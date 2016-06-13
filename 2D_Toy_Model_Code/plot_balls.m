%% Spectral Clustering
ball_clustering = load('ball_clustering_8176.txt');
m = size(ball_clustering, 1);
a = zeros(m, 3);
figure; hold on
for i=1:m
  a(i, 1) = ball_clustering(i, 3);
  a(i, 2) = log(abs(ball_clustering(i, 4)));
  a(i, 3) = ball_clustering(i, 5)/abs(ball_clustering(i, 4));
  plot(ball_clustering(i, 1), ball_clustering(i, 2),'ro', 'MarkerSize', 1.0);
end
s = 0.2; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(ball_clustering(:, 1), ball_clustering(:, 2), 1, a(:, 3), 'Linewidth', 1.0); colormap(jet), colorbar
%axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x');
ylabel('y');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Regular Plot
total_weight = load('total_weight_on_each_ball_8176.txt');
m = size(total_weight, 1);
a = zeros(m, 3);
figure; hold on
for i=1:m
    if total_weight(i, 3) > 0.0
        a(i, 1) = total_weight(i, 1);
        a(i, 2) = total_weight(i, 2);
        a(i, 3) = total_weight(i, 3);
        plot(total_weight(i, 1), total_weight(i, 2), 'ro', 'MarkerSize', 1.0);
    end
end
s = 0.2; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(a(:, 1), a(:, 2), log(a(:, 3)), 'Linewidth', 1.0); colormap(jet), colorbar, caxis([-20 0])
%axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x')
ylabel('y')
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Initial Conditions
x = (-0.6);
y = (0.0);
a = (0.0);
figure; hold on
plot(x, y, 'ro', 'MarkerSize', 1.0);
s = 0.2; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(x, y, 1, a, 'Linewidth', 1.0); colormap(jet), colorbar
axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x')
ylabel('y')
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);
