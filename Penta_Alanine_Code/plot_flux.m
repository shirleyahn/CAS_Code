%%
load('flux_total_40.txt')
load('flux_total_60.txt')
load('flux_total_60_sc.txt')
load('flux_total_60_sc2.txt')
load('flux_total_80.txt')
load('flux_total_100.txt')
load('flux_total_120.txt')

figure;hold on;
plot((1:300)*0.5, ones(300,1)*0.033, '-.k');
plot(flux_total_40(:,1),flux_total_40(:,2),'r');
plot(flux_total_60(:,1),flux_total_60(:,2),'y');
plot(flux_total_60_sc(:,1),flux_total_60_sc(:,2),'g');
plot(flux_total_60_sc2(:,1),flux_total_60_sc2(:,2),'b');
plot(flux_total_80(:,1),flux_total_80(:,2),'m');
plot(flux_total_100(:,1),flux_total_100(:,2),'c');
plot(flux_total_120(:,1),flux_total_120(:,2),'k');
legend('exact','40','60','60 sc','60 sc2','80','100','120');


%%
figure;hold on;
plot((1:300)*0.5, ones(300,1)*0.054, '-.k');
plot(flux_total_40(:,1),flux_total_40(:,3),'r');
plot(flux_total_60(:,1),flux_total_60(:,3),'y');
plot(flux_total_60_sc(:,1),flux_total_60_sc(:,3),'g');
plot(flux_total_60_sc2(:,1),flux_total_60_sc2(:,3),'b');
plot(flux_total_80(:,1),flux_total_80(:,3),'m');
plot(flux_total_100(:,1),flux_total_100(:,3),'c');
plot(flux_total_120(:,1),flux_total_120(:,3),'k');
legend('exact','40','60','60 sc','60 sc2','80','100','120');
