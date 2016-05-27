%%
load('flux_total_80.txt')
load('flux_total_80_sc.txt')

figure;hold on;
plot(flux_total_80(:,1), ones(153,1)*0.033, 'k');
plot(flux_total_80(:,1),flux_total_80(:,2),'r');
plot(flux_total_80_sc(:,1),flux_total_80_sc(:,2),'g');


%%
figure;hold on;
plot(flux_total_80(:,1), ones(153,1)*0.054, 'k');
plot(flux_total_80(:,1),flux_total_80(:,3),'r');
plot(flux_total_80_sc(:,1),flux_total_80_sc(:,3),'g');
