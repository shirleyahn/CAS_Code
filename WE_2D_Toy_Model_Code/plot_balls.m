load('total_weight_on_each_ball_4.txt')
total_weight=total_weight_on_each_ball_4;
m=size(total_weight,1);
a=zeros(m,1);
c=zeros(m,3);
figure; hold on
for i=1:m
    a(i) = total_weight(i,3);
    plot(total_weight(i,1),total_weight(i,2),'ro', 'MarkerSize', 1.0);
end
s = 0.2; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(total_weight(:,1),total_weight(:,2),1,log(a),'Linewidth',1.0),
colormap(jet),colorbar,caxis([-13 -3]);
axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x')
ylabel('y')
% Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);
