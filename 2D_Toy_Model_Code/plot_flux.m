load('flux_avg_1ball_10steps_188.txt')
%load('flux_avg_1ball_10steps_193.txt')
%load('flux_avg_1ball_10steps_200.txt')
%load('flux_avg_1ball_10steps_215.txt')
%load('flux_avg_1ball_10steps_220.txt')
load('flux_avg_1ball_10steps_sc_200_2000.txt')
%load('flux_avg_1ball_10steps_sc_200_4000.txt')
%load('flux_avg_1ball_10steps_sc_200_4000_no_od.txt')
%load('flux_avg_1ball_10steps_sc_250_2500.txt')
%load('flux_avg_1ball_10steps_sc_250_5000.txt')
flux_avg_1ball_10steps_r_0_05 = load('flux_total_1ball_10steps_r_0.05_1.txt');
flux_avg_1ball_10steps_r_0_1 = load('flux_avg_1ball_10steps_r_0.1.txt');
flux_avg_1ball_10steps_r_0_2 = load('flux_avg_1ball_10steps_r_0.2.txt');
flux_avg_1ball_10steps_r_0_4 = load('flux_avg_1ball_10steps_r_0.4.txt');
%load('flux_total_1ball_10steps_sc_200_10_1.txt')
%load('flux_total_1ball_10steps_sc_200_10_no_od_1.txt')
%load('flux_total_1ball_10steps_sc_200_20_1.txt')
%load('flux_total_1ball_10steps_sc_200_20_no_od_1.txt')
%load('flux_total_1ball_10steps_sc_200_40_1.txt')
%load('flux_total_1ball_10steps_sc_200_40_no_od_1.txt')
%load('flux_total_1ball_10steps_sc_200_50_1.txt')
%load('flux_total_1ball_10steps_sc_200_50_no_od_1.txt')
%load('flux_total_1ball_10steps_sc_200_100_1.txt')
%load('flux_total_1ball_10steps_sc_200_200_1.txt')
%load('flux_total_1ball_10steps_sc_200_400_1.txt')

load('flux_std_1ball_10steps_188.txt')
%load('flux_std_1ball_10steps_193.txt')
%load('flux_std_1ball_10steps_200.txt')
%load('flux_std_1ball_10steps_215.txt')
%load('flux_std_1ball_10steps_220.txt')
load('flux_std_1ball_10steps_sc_200_2000.txt')
%load('flux_std_1ball_10steps_sc_200_4000.txt')
%load('flux_std_1ball_10steps_sc_200_4000_no_od.txt')
%load('flux_std_1ball_10steps_sc_250_2500.txt')
%load('flux_std_1ball_10steps_sc_250_5000.txt')
flux_std_1ball_10steps_r_0_1 = load('flux_std_1ball_10steps_r_0.1.txt');
flux_std_1ball_10steps_r_0_2 = load('flux_std_1ball_10steps_r_0.2.txt');
flux_std_1ball_10steps_r_0_4 = load('flux_std_1ball_10steps_r_0.4.txt');

total_weight_r_0_05 = load('CAS_1ball_10steps_r_0.05_1/total_weight.txt');
total_weight_r_0_1 = load('total_weight_avg_CAS_1ball_10steps_r_0.1.txt');
total_weight_r_0_2 = load('total_weight_avg_CAS_1ball_10steps_r_0.2.txt');
total_weight_r_0_4 = load('total_weight_avg_CAS_1ball_10steps_r_0.4.txt');
total_weight_188 = load('total_weight_avg_CAS_1ball_10steps_188.txt');

figure;hold on;
plot((1:1:5000),ones(5000,1)*1.1691698347e-04,'-.k');
plot(flux_avg_1ball_10steps_r_0_05(:,1),flux_avg_1ball_10steps_r_0_05(:,2),'r');
plot(flux_avg_1ball_10steps_r_0_4(:,1),flux_avg_1ball_10steps_r_0_4(:,2),'g');
plot(flux_avg_1ball_10steps_r_0_1(:,1),flux_avg_1ball_10steps_r_0_1(:,2),'b');
plot(flux_avg_1ball_10steps_188(:,1),flux_avg_1ball_10steps_188(:,2),'c');
plot(flux_avg_1ball_10steps_sc_200_2000(:,1),flux_avg_1ball_10steps_sc_200_2000(:,2),'m');

xlabel('time (# of steps)')
ylabel('forward flux')
axis([1, 5000, 0, 2.5e-04])
legend('exact','r=0.05 unbounded (avg: 634 balls)', 'r=0.4 unbounded (avg: 14 balls)', 'r=0.1 unbounded (avg: 251 balls)', 'r=0.1 bounded (avg: 188 balls)','spectral clustering with 200 balls, 2000 walkers (avg: 188 balls)')
%%
figure;hold on;
plot((1:1:5000),ones(5000,1)*1.1691698347e-04,'-.k');
plot(flux_avg_1ball_10steps_r_0_05(:,1),flux_avg_1ball_10steps_r_0_05(:,3),'r');
plot(flux_avg_1ball_10steps_r_0_4(:,1),flux_avg_1ball_10steps_r_0_4(:,3),'g');
plot(flux_avg_1ball_10steps_r_0_1(:,1),flux_avg_1ball_10steps_r_0_1(:,3),'b');
plot(flux_avg_1ball_10steps_188(:,1),flux_avg_1ball_10steps_188(:,3),'c');
plot(flux_avg_1ball_10steps_sc_200_2000(:,1),flux_avg_1ball_10steps_sc_200_2000(:,3),'m');

xlabel('time (# of steps)')
ylabel('backward flux')
axis([1, 5000, 0, 2.5e-04])
legend('exact','r=0.05 unbounded (avg: 634 balls)', 'r=0.4 unbounded (avg: 14 balls)', 'r=0.1 unbounded (avg: 251 balls)', 'r=0.1 bounded (avg: 188 balls)','spectral clustering with 200 balls, 2000 walkers (avg: 188 balls)')
