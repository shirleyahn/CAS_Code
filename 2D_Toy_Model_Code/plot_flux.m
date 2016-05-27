load('flux_avg_1ball_10steps_188.txt')
load('flux_avg_1ball_10steps_193.txt')
load('flux_avg_1ball_10steps_200.txt')
load('flux_avg_1ball_10steps_215.txt')
load('flux_avg_1ball_10steps_220.txt')
load('flux_avg_1ball_10steps_sc_200_400.txt')
load('flux_avg_1ball_10steps_sc_200_500.txt')
load('flux_avg_1ball_10steps_sc_200_800.txt')
load('flux_avg_1ball_10steps_sc_200_1000.txt')
load('flux_avg_1ball_10steps_sc_200_1500.txt')
load('flux_avg_1ball_10steps_sc_200_2000.txt')
load('flux_avg_1ball_10steps_sc_200_4000.txt')
load('flux_avg_1ball_10steps_sc_200_4000_no_od.txt')
load('flux_avg_1ball_10steps_sc_250_2500.txt')
load('flux_avg_1ball_10steps_sc_250_5000.txt')
flux_avg_1ball_10steps_r_0_05 = load('flux_avg_1ball_10steps_r_0.05.txt');
flux_avg_1ball_10steps_r_0_06 = load('flux_avg_1ball_10steps_r_0.06.txt');
flux_avg_1ball_10steps_r_0_08 = load('flux_avg_1ball_10steps_r_0.08.txt');
flux_avg_1ball_10steps_r_0_1 = load('flux_avg_1ball_10steps_r_0.1.txt');
flux_avg_1ball_10steps_r_0_2 = load('flux_avg_1ball_10steps_r_0.2.txt');
flux_avg_1ball_10steps_r_0_4 = load('flux_avg_1ball_10steps_r_0.4.txt');
flux_avg_1ball_1step_r_0_1_beta_10 = load('flux_avg_1ball_1step_r_0.1_beta_10.txt');
flux_avg_1ball_1step_r_0_1_beta_20 = load('flux_avg_1ball_1step_r_0.1_beta_20.txt');
flux_avg_1ball_1step_r_0_1_beta_30 = load('flux_avg_1ball_1step_r_0.1_beta_30.txt');
flux_avg_1ball_1step_r_0_2_beta_10 = load('flux_avg_1ball_1step_r_0.2_beta_10.txt');
flux_avg_1ball_1step_r_0_2_beta_20 = load('flux_avg_1ball_1step_r_0.2_beta_20.txt');
flux_avg_1ball_1step_r_0_2_beta_30 = load('flux_avg_1ball_1step_r_0.2_beta_30.txt');
flux_avg_1ball_1step_r_0_4_beta_10 = load('flux_avg_1ball_1step_r_0.4_beta_10.txt');
flux_avg_1ball_1step_r_0_4_beta_20 = load('flux_avg_1ball_1step_r_0.4_beta_20.txt');
flux_avg_1ball_1step_r_0_4_beta_30 = load('flux_avg_1ball_1step_r_0.4_beta_30.txt');

load('flux_std_1ball_10steps_188.txt')
load('flux_std_1ball_10steps_193.txt')
load('flux_std_1ball_10steps_200.txt')
load('flux_std_1ball_10steps_215.txt')
load('flux_std_1ball_10steps_220.txt')
load('flux_std_1ball_10steps_sc_200_400.txt')
load('flux_std_1ball_10steps_sc_200_500.txt')
load('flux_std_1ball_10steps_sc_200_800.txt')
load('flux_std_1ball_10steps_sc_200_1000.txt')
load('flux_std_1ball_10steps_sc_200_1500.txt')
load('flux_std_1ball_10steps_sc_200_2000.txt')
load('flux_std_1ball_10steps_sc_200_4000.txt')
load('flux_std_1ball_10steps_sc_200_4000_no_od.txt')
load('flux_std_1ball_10steps_sc_250_2500.txt')
load('flux_std_1ball_10steps_sc_250_5000.txt')
flux_std_1ball_10steps_r_0_05 = load('flux_std_1ball_10steps_r_0.05.txt');
flux_std_1ball_10steps_r_0_06 = load('flux_std_1ball_10steps_r_0.06.txt');
flux_std_1ball_10steps_r_0_08 = load('flux_std_1ball_10steps_r_0.08.txt');
flux_std_1ball_10steps_r_0_1 = load('flux_std_1ball_10steps_r_0.1.txt');
flux_std_1ball_10steps_r_0_2 = load('flux_std_1ball_10steps_r_0.2.txt');
flux_std_1ball_10steps_r_0_4 = load('flux_std_1ball_10steps_r_0.4.txt');
flux_std_1ball_1step_r_0_1_beta_10 = load('flux_std_1ball_1step_r_0.1_beta_10.txt');
flux_std_1ball_1step_r_0_1_beta_20 = load('flux_std_1ball_1step_r_0.1_beta_20.txt');
flux_std_1ball_1step_r_0_1_beta_30 = load('flux_std_1ball_1step_r_0.1_beta_30.txt');
flux_std_1ball_1step_r_0_2_beta_10 = load('flux_std_1ball_1step_r_0.2_beta_10.txt');
flux_std_1ball_1step_r_0_2_beta_20 = load('flux_std_1ball_1step_r_0.2_beta_20.txt');
flux_std_1ball_1step_r_0_2_beta_30 = load('flux_std_1ball_1step_r_0.2_beta_30.txt');
flux_std_1ball_1step_r_0_4_beta_10 = load('flux_std_1ball_1step_r_0.4_beta_10.txt');
flux_std_1ball_1step_r_0_4_beta_20 = load('flux_std_1ball_1step_r_0.4_beta_20.txt');
flux_std_1ball_1step_r_0_4_beta_30 = load('flux_std_1ball_1step_r_0.4_beta_30.txt');


%%
figure;hold on;
start = 1;
endpt = 500;
freq = 1;
plot((start:1:endpt),ones((endpt-start)+1,1)*1.1691698347e-05,'-.k');
%plot((start:1:endpt),ones((endpt-start)+1,1)*4.8006579851e-08,'-.k');
%plot((start:1:endpt),ones((endpt-start)+1,1)*1.3137196888e-10,'-.k');
%errorbar(flux_avg_1ball_10steps_r_0_05(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_05(start:freq:endpt,2),flux_std_1ball_10steps_r_0_05(start:freq:endpt,2),'r');
%errorbar(flux_avg_1ball_10steps_r_0_06(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_06(start:freq:endpt,2),flux_std_1ball_10steps_r_0_06(start:freq:endpt,2),'y');
%errorbar(flux_avg_1ball_10steps_r_0_08(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_08(start:freq:endpt,2),flux_std_1ball_10steps_r_0_08(start:freq:endpt,2),'g');
%errorbar(flux_avg_1ball_10steps_r_0_1(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_1(start:freq:endpt,2),flux_std_1ball_10steps_r_0_1(start:freq:endpt,2),'b');
%errorbar(flux_avg_1ball_10steps_r_0_2(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_2(start:freq:endpt,2),flux_std_1ball_10steps_r_0_2(start:freq:endpt,2),'m');
%errorbar(flux_avg_1ball_10steps_r_0_4(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_4(start:freq:endpt,2),flux_std_1ball_10steps_r_0_4(start:freq:endpt,2),'c');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_10(start:freq:endpt,2),flux_std_1ball_1step_r_0_1_beta_10(start:freq:endpt,2),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_10(start:freq:endpt,2),flux_std_1ball_1step_r_0_2_beta_10(start:freq:endpt,2),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_10(start:freq:endpt,2),flux_std_1ball_1step_r_0_4_beta_10(start:freq:endpt,2),'b');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_20(start:freq:endpt,2),flux_std_1ball_1step_r_0_1_beta_20(start:freq:endpt,2),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_20(start:freq:endpt,2),flux_std_1ball_1step_r_0_2_beta_20(start:freq:endpt,2),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_20(start:freq:endpt,2),flux_std_1ball_1step_r_0_4_beta_20(start:freq:endpt,2),'b');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_30(start:freq:endpt,2),flux_std_1ball_1step_r_0_1_beta_30(start:freq:endpt,2),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_30(start:freq:endpt,2),flux_std_1ball_1step_r_0_2_beta_30(start:freq:endpt,2),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_30(start:freq:endpt,2),flux_std_1ball_1step_r_0_4_beta_30(start:freq:endpt,2),'b');
%errorbar(flux_avg_1ball_10steps_sc_200_400(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_400(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_400(start:freq:endpt,2),'r');
%errorbar(flux_avg_1ball_10steps_sc_200_500(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_500(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_500(start:freq:endpt,2),'y');
%errorbar(flux_avg_1ball_10steps_sc_200_800(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_800(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_800(start:freq:endpt,2),'g');
%errorbar(flux_avg_1ball_10steps_sc_200_1000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_1000(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_1000(start:freq:endpt,2),'b');
%errorbar(flux_avg_1ball_10steps_sc_200_1500(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_1500(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_1500(start:freq:endpt,2),'m');
%errorbar(flux_avg_1ball_10steps_sc_200_2000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_2000(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_2000(start:freq:endpt,2),'c');
%errorbar(flux_avg_1ball_10steps_sc_200_4000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_4000(start:freq:endpt,2),flux_std_1ball_10steps_sc_200_4000(start:freq:endpt,2),'k');

xlabel('time (# of steps)')
ylabel('forward flux')
%axis([start, endpt, 1.16e-04, 1.2e-04])
%legend('exact','r=0.08 unbounded (avg: 334 balls)','r=0.1 unbounded (avg: 251 balls)', 'r=0.1 bounded (avg: 193 balls)','spectral clustering with 200 balls, 4000 walkers (avg: 193 balls)')


%%
figure;hold on;
start = 1;
endpt = 500;
freq = 1;
plot((start:1:endpt),ones((endpt-start)+1,1)*1.1691698347e-05,'-.k');
%plot((start:1:endpt),ones((endpt-start)+1,1)*4.8006579851e-08,'-.k');
%plot((start:1:endpt),ones((endpt-start)+1,1)*1.3137196888e-10,'-.k');
%errorbar(flux_avg_1ball_10steps_r_0_05(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_05(start:freq:endpt,3),flux_std_1ball_10steps_r_0_05(start:freq:endpt,3),'r');
%errorbar(flux_avg_1ball_10steps_r_0_06(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_06(start:freq:endpt,3),flux_std_1ball_10steps_r_0_06(start:freq:endpt,3),'y');
%errorbar(flux_avg_1ball_10steps_r_0_08(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_08(start:freq:endpt,3),flux_std_1ball_10steps_r_0_08(start:freq:endpt,3),'g');
%errorbar(flux_avg_1ball_10steps_r_0_1(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_1(start:freq:endpt,3),flux_std_1ball_10steps_r_0_1(start:freq:endpt,3),'b');
%errorbar(flux_avg_1ball_10steps_r_0_2(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_2(start:freq:endpt,3),flux_std_1ball_10steps_r_0_2(start:freq:endpt,3),'m');
%errorbar(flux_avg_1ball_10steps_r_0_4(start:freq:endpt,1),flux_avg_1ball_10steps_r_0_4(start:freq:endpt,3),flux_std_1ball_10steps_r_0_4(start:freq:endpt,3),'c');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_10(start:freq:endpt,3),flux_std_1ball_1step_r_0_1_beta_10(start:freq:endpt,3),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_10(start:freq:endpt,3),flux_std_1ball_1step_r_0_2_beta_10(start:freq:endpt,3),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_10(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_10(start:freq:endpt,3),flux_std_1ball_1step_r_0_4_beta_10(start:freq:endpt,3),'b');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_20(start:freq:endpt,3),flux_std_1ball_1step_r_0_1_beta_20(start:freq:endpt,3),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_20(start:freq:endpt,3),flux_std_1ball_1step_r_0_2_beta_20(start:freq:endpt,3),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_20(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_20(start:freq:endpt,3),flux_std_1ball_1step_r_0_4_beta_20(start:freq:endpt,3),'b');
%errorbar(flux_avg_1ball_1step_r_0_1_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_1_beta_30(start:freq:endpt,3),flux_std_1ball_1step_r_0_1_beta_30(start:freq:endpt,3),'r');
%errorbar(flux_avg_1ball_1step_r_0_2_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_2_beta_30(start:freq:endpt,3),flux_std_1ball_1step_r_0_2_beta_30(start:freq:endpt,3),'g');
%errorbar(flux_avg_1ball_1step_r_0_4_beta_30(start:freq:endpt,1),flux_avg_1ball_1step_r_0_4_beta_30(start:freq:endpt,3),flux_std_1ball_1step_r_0_4_beta_30(start:freq:endpt,3),'b');
%errorbar(flux_avg_1ball_10steps_sc_200_400(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_400(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_400(start:freq:endpt,3),'r');
%errorbar(flux_avg_1ball_10steps_sc_200_500(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_500(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_500(start:freq:endpt,3),'y');
%errorbar(flux_avg_1ball_10steps_sc_200_800(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_800(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_800(start:freq:endpt,3),'g');
%errorbar(flux_avg_1ball_10steps_sc_200_1000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_1000(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_1000(start:freq:endpt,3),'b');
%errorbar(flux_avg_1ball_10steps_sc_200_1500(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_1500(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_1500(start:freq:endpt,3),'m');
%errorbar(flux_avg_1ball_10steps_sc_200_2000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_2000(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_2000(start:freq:endpt,3),'c');
%errorbar(flux_avg_1ball_10steps_sc_200_4000(start:freq:endpt,1),flux_avg_1ball_10steps_sc_200_4000(start:freq:endpt,3),flux_std_1ball_10steps_sc_200_4000(start:freq:endpt,3),'k');

xlabel('time (# of steps)')
ylabel('backward flux')
%axis([4500, 5000, 1.16e-04, 1.2e-04])
%legend('exact','r=0.08 unbounded (avg: 334 balls)','r=0.1 unbounded (avg: 251 balls)', 'r=0.1 bounded (avg: 193 balls)','spectral clustering with 200 balls, 4000 walkers (avg: 193 balls)')


%%
load('flux_avg_1ball_1step_beta_20_sc_200_2000_20.txt')
load('flux_avg_1ball_1step_beta_20_sc_200_4000_20.txt')
flux_avg_1ball_1step_beta_20_r_0_1 = load('flux_avg_1ball_1step_beta_20_r_0.1.txt');
flux_avg_1ball_1step_beta_20_r_0_2 = load('flux_avg_1ball_1step_beta_20_r_0.2.txt');
flux_avg_1ball_1step_beta_20_r_0_4 = load('flux_avg_1ball_1step_beta_20_r_0.4.txt');
flux_avg_1ball_1step_beta_20_r_0_8 = load('flux_avg_1ball_1step_beta_20_r_0.8.txt');

load('flux_std_1ball_1step_beta_20_sc_200_2000_20.txt')
load('flux_std_1ball_1step_beta_20_sc_200_4000_20.txt')
flux_std_1ball_1step_beta_20_r_0_1 = load('flux_std_1ball_1step_beta_20_r_0.1.txt');
flux_std_1ball_1step_beta_20_r_0_2 = load('flux_std_1ball_1step_beta_20_r_0.2.txt');
flux_std_1ball_1step_beta_20_r_0_4 = load('flux_std_1ball_1step_beta_20_r_0.4.txt');
flux_std_1ball_1step_beta_20_r_0_8 = load('flux_std_1ball_1step_beta_20_r_0.8.txt');

load('total_time_avg_1ball_1step_beta_20_r_0.1.txt')
load('total_time_avg_1ball_1step_beta_20_r_0.2.txt')
load('total_time_avg_1ball_1step_beta_20_r_0.4.txt')
load('total_time_avg_1ball_1step_beta_20_r_0.8.txt')
load('total_time_avg_1ball_1step_beta_20_sc_200_2000_20.txt')


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*4.8006567250e-08,'-.k');
%plot((1:4000),ones(4000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_20_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_1(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_1(1:freq:end,2),'r');
errorbar(flux_avg_1ball_1step_beta_20_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_2(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_2(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_4(1:freq:end,2),'b');
errorbar(flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_8(1:freq:end,2),'m');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,2),'b');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.8 unbounded','r = 0.1 spectral clustering')


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*4.8006567250e-08,'-.k');
%plot((1:4000),ones(4000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_20_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_1(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_1(1:freq:end,3),'r');
errorbar(flux_avg_1ball_1step_beta_20_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_2(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_2(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_4(1:freq:end,3),'b');
errorbar(flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_8(1:freq:end,3),'m');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,3),'b');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('backward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.8 unbounded','r = 0.1 spectral clustering')


%%
flux_avg_1ball_1step_beta_4_r_0_1 = load('flux_avg_1ball_1step_beta_4_r_0.1.txt');
flux_avg_1ball_1step_beta_4_r_0_2 = load('flux_avg_1ball_1step_beta_4_r_0.2.txt');
flux_avg_1ball_1step_beta_4_r_0_4 = load('flux_avg_1ball_1step_beta_4_r_0.4.txt');

flux_std_1ball_1step_beta_4_r_0_1 = load('flux_std_1ball_1step_beta_4_r_0.1.txt');
flux_std_1ball_1step_beta_4_r_0_2 = load('flux_std_1ball_1step_beta_4_r_0.2.txt');
flux_std_1ball_1step_beta_4_r_0_4 = load('flux_std_1ball_1step_beta_4_r_0.4.txt');

load('total_time_avg_1ball_1step_beta_4_r_0.1.txt')
load('total_time_avg_1ball_1step_beta_4_r_0.2.txt')
load('total_time_avg_1ball_1step_beta_4_r_0.4.txt')

%%
figure; hold on;
freq = 100;
%plot((1:4000),ones(4000,1)*4.8006567250e-08,'-.k');
plot((1:4000),ones(4000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_1(1:freq:end,2),'r');
errorbar(flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_2(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_4(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,2),flux_std_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,2),'m');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')


%%
figure; hold on;
freq = 100;
%plot((1:4000),ones(4000,1)*4.8006567250e-08,'-.k');
plot((1:4000),ones(4000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_1(1:freq:end,3),'r');
errorbar(flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_2(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_4(1:freq:end,3),'b');
%errorbar(flux_avg_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,3),flux_std_1ball_1step_beta_4_sc_400_8000_20(1:freq:end,3),'m');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')
