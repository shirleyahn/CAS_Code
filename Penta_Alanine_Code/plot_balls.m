%% Penta Alanine Spectral Clustering
ball_clustering = load('ball_clustering_16.txt');
m = size(ball_clustering, 1);
a = zeros(m, 3);
figure; hold on
for i=1:m
  a(i, 1) = ball_clustering(i, 7);
  a(i, 2) = abs(ball_clustering(i, 8));
  a(i, 3) = ball_clustering(i, 9)/abs(ball_clustering(i, 8));
  plot(ball_clustering(i, 5), ball_clustering(i, 6),'ro', 'MarkerSize', 1.0);
end
s = 30.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(ball_clustering(:, 5), ball_clustering(:, 6), 1, a(:, 1), 'Linewidth', 1.0), colormap(jet), colorbar;
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
total_weight = load('total_weight_on_each_ball_161.txt');
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 1), total_weight(i, 2), 'ro', 'MarkerSize', 1.0);
end
s = 20.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 1), total_weight(:, 2), 1, -0.0019872041*300*log(a), 'Linewidth', 10.0), colormap(jet), colorbar;
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
total_weight = load('total_weight_on_each_ball_2.txt');
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 3), total_weight(i, 4), 'ro', 'MarkerSize', 1.0);
end
s = 30.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 3), total_weight(:, 4), 1, -0.0019872041*300*log(a), 'Linewidth', 10.0), colormap(jet), colorbar;
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
total_weight = load('total_weight_on_each_ball_2.txt');
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
    a(i) = total_weight(i, 7);
    plot(total_weight(i, 5), total_weight(i, 6), 'ro', 'MarkerSize', 1.0);
end
s = 30.0; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 5), total_weight(:, 6), 1, -0.0019872041*300*log(a), 'Linewidth', 10.0), colormap(jet), colorbar;
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


%% Penta Alanine Graph
A = load('transition_matrix_37.txt');
B = load('balls_37.txt');
initial = [];
final = [];
for i = 1:size(B, 1)
    if -100.0 <= B(i, 1) && B(i, 1) <= -30.0 && -90.0 <= B(i, 2) && ... 
            B(i, 2) <= -10.0 && -100.0 <= B(i, 3) && B(i, 3) <= -30.0 ...
            && -90.0 <= B(i, 4) && B(i, 4) <= -10.0 && -100.0 <= B(i, 5) ...
            && B(i, 5) <= -30.0  && -90.0 <= B(i, 6) && B(i, 6) <= -10.0
        initial = [initial i];
    elseif -180.0 <= B(i, 1) && B(i, 1) <= -55.0 && ((105.0 <= B(i, 2) && ...
            B(i, 2) <= 180.0) || (-180.0 <= B(i, 2) && B(i, 2) <= -155.0)) ...
            && -180.0 <= B(i, 3) && B(i, 3) <= -55.0 && ((105.0 <= B(i, 4) && ...
            B(i, 4) <= 180.0) || (-180.0 <= B(i, 4) && B(i, 4) <= -155.0)) ...
            && -180.0 <= B(i, 5) && B(i, 5) <= -55.0 && ((105.0 <= B(i, 6) ...
            && B(i, 6) <= 180.0) || (-180.0 <= B(i, 6) && B(i, 6) <= -155.0))
        final = [final i];
    end
end
%G = digraph(A, 'OmitSelfLoops');
%LWidths = 5*-log(G.Edges.Weight)/max(-log(G.Edges.Weight));
k = 0;
figure; hold on;
%for i = 1:size(A, 1)
    %for j = 1:size(A, 2)
        %if i ~= j && A(i, j) ~= 0.0
            %k = k + 1;
            %if LWidths(k) > 3
                %p1 = [B(i, 1), B(i, 2)];
                %p2 = [B(j, 1), B(j, 2)];
                %dp = p2-p1;
                %quiver(p1(1), p1(2), dp(1), dp(2), 'LineWidth', LWidths(k))
                %text(p1(1), p1(2), num2str(log(A(i,j))))
            %end
        %end 
    %end
%end
for i = 1:size(initial, 2)
    for j = 1:size(final, 2)
        ii = initial(i);
        jj = final(j);
        p1 = [B(ii, 1), B(ii, 2)];
        p2 = [B(jj, 1), B(jj, 2)];
        dp = p2-p1;
        quiver(p1(1), p1(2), dp(1), dp(2))
        text(p2(1), p2(2), num2str(A(ii,jj)))
    end
end    
axis([-180 180 -180 180]);