load('flux_avg_1ball_1step_beta_4_sc_600_12000_20.txt')
flux_avg_1ball_1step_beta_4_r_0_1 = load('flux_avg_1ball_1step_beta_4_r_0.1.txt');
flux_avg_1ball_1step_beta_4_r_0_2 = load('flux_avg_1ball_1step_beta_4_r_0.2.txt');
flux_avg_1ball_1step_beta_4_r_0_4 = load('flux_avg_1ball_1step_beta_4_r_0.4.txt');

load('flux_std_1ball_1step_beta_4_sc_600_12000_20.txt')
flux_std_1ball_1step_beta_4_r_0_1 = load('flux_std_1ball_1step_beta_4_r_0.1.txt');
flux_std_1ball_1step_beta_4_r_0_2 = load('flux_std_1ball_1step_beta_4_r_0.2.txt');
flux_std_1ball_1step_beta_4_r_0_4 = load('flux_std_1ball_1step_beta_4_r_0.4.txt');


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_1(1:freq:end,2),'r');
errorbar(flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_2(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,2),flux_std_1ball_1step_beta_4_r_0_4(1:freq:end,2),'b');
errorbar(flux_avg_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,2),flux_std_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,2),'m');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*2.4603069588e-10,'-.k');
errorbar(flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_1(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_1(1:freq:end,3),'r');
errorbar(flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_2(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_2(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_4_r_0_4(1:freq:end,3),flux_std_1ball_1step_beta_4_r_0_4(1:freq:end,3),'b');
errorbar(flux_avg_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,3),flux_std_1ball_1step_beta_4_sc_600_12000_20(1:freq:end,3),'m');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('backward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')


%%
load('flux_avg_1ball_1step_beta_15_sc_200_4000_20.txt')
load('flux_avg_1ball_1step_beta_15_sc_200_4000_20_v2.txt')
load('flux_avg_2balls_1step_beta_15_sc_200_4000_20.txt')
load('flux_avg_1ball_1step_beta_15_sc_400_8000_20.txt')
load('flux_avg_1ball_1step_beta_15_sc_600_12000_20.txt')
load('flux_avg_1ball_1step_beta_15_sc_800_16000_20.txt')
load('flux_avg_1ball_1step_beta_15_sc_1000_20000_20.txt')
flux_avg_1ball_1step_beta_15_r_0_1 = load('flux_avg_1ball_1step_beta_15_r_0.1.txt');
flux_avg_1ball_1step_beta_15_r_0_8 = load('flux_avg_1ball_1step_beta_15_r_0.8.txt');
flux_avg_1ball_1step_beta_15_r_0_8_v2 = load('flux_avg_1ball_1step_beta_15_r_0.8_v2.txt');
flux_avg_2balls_1step_beta_15_r_0_1 = load('flux_avg_2balls_1step_beta_15_r_0.1.txt');
flux_avg_2balls_1step_beta_15_r_0_8 = load('flux_avg_2balls_1step_beta_15_r_0.8.txt');

load('flux_std_1ball_1step_beta_15_sc_200_4000_20.txt')
load('flux_std_1ball_1step_beta_15_sc_200_4000_20_v2.txt')
load('flux_std_2balls_1step_beta_15_sc_200_4000_20.txt')
load('flux_std_1ball_1step_beta_15_sc_400_8000_20.txt')
load('flux_std_1ball_1step_beta_15_sc_600_12000_20.txt')
load('flux_std_1ball_1step_beta_15_sc_800_16000_20.txt')
load('flux_std_1ball_1step_beta_15_sc_1000_20000_20.txt')
flux_std_1ball_1step_beta_15_r_0_1 = load('flux_std_1ball_1step_beta_15_r_0.1.txt');
flux_std_1ball_1step_beta_15_r_0_8 = load('flux_std_1ball_1step_beta_15_r_0.8.txt');
flux_std_1ball_1step_beta_15_r_0_8_v2 = load('flux_std_1ball_1step_beta_15_r_0.8_v2.txt');
flux_std_2balls_1step_beta_15_r_0_1 = load('flux_std_2balls_1step_beta_15_r_0.1.txt');
flux_std_2balls_1step_beta_15_r_0_8 = load('flux_std_2balls_1step_beta_15_r_0.8.txt');


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*8.1443713035e-07,'-.k');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_1(1:freq:end,2),flux_std_1ball_1step_beta_15_r_0_1(1:freq:end,2),'r');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_8(1:freq:end,2),flux_std_1ball_1step_beta_15_r_0_8(1:freq:end,2),'g');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_8_v2(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_8_v2(1:freq:end,2),flux_std_1ball_1step_beta_15_r_0_8_v2(1:freq:end,2),'b');
%errorbar(flux_avg_2balls_1step_beta_15_r_0_1(1:freq:end,1),flux_avg_2balls_1step_beta_15_r_0_1(1:freq:end,2),flux_std_2balls_1step_beta_15_r_0_1(1:freq:end,2),'r');
%errorbar(flux_avg_2balls_1step_beta_15_r_0_8(1:freq:end,1),flux_avg_2balls_1step_beta_15_r_0_8(1:freq:end,2),flux_std_2balls_1step_beta_15_r_0_8(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,2),'r');
errorbar(flux_avg_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,2),'y');
errorbar(flux_avg_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,1),flux_avg_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,2),flux_std_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,2),'b');
errorbar(flux_avg_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,2),'m');
errorbar(flux_avg_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,2),'c');
errorbar(flux_avg_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,2),flux_std_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,2),'k');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')


%%
figure; hold on;
freq = 100;
plot((1:20000),ones(20000,1)*8.1443713035e-07,'-.k');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_1(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_1(1:freq:end,2),flux_std_1ball_1step_beta_15_r_0_1(1:freq:end,2),'r');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_8(1:freq:end,3),flux_std_1ball_1step_beta_15_r_0_8(1:freq:end,3),'g');
%errorbar(flux_avg_1ball_1step_beta_15_r_0_8_v2(1:freq:end,1),flux_avg_1ball_1step_beta_15_r_0_8_v2(1:freq:end,3),flux_std_1ball_1step_beta_15_r_0_8_v2(1:freq:end,3),'b');
%errorbar(flux_avg_2balls_1step_beta_15_r_0_1(1:freq:end,1),flux_avg_2balls_1step_beta_15_r_0_1(1:freq:end,3),flux_std_2balls_1step_beta_15_r_0_1(1:freq:end,3),'r');
%errorbar(flux_avg_2balls_1step_beta_15_r_0_8(1:freq:end,1),flux_avg_2balls_1step_beta_15_r_0_8(1:freq:end,3),flux_std_2balls_1step_beta_15_r_0_8(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_200_4000_20(1:freq:end,3),'r');
errorbar(flux_avg_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_200_4000_20_v2(1:freq:end,3),'y');
errorbar(flux_avg_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,1),flux_avg_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,3),flux_std_2balls_1step_beta_15_sc_200_4000_20(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_400_8000_20(1:freq:end,3),'b');
errorbar(flux_avg_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_600_12000_20(1:freq:end,3),'m');
errorbar(flux_avg_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_800_16000_20(1:freq:end,3),'c');
errorbar(flux_avg_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,1),flux_avg_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,3),flux_std_1ball_1step_beta_15_sc_1000_20000_20(1:freq:end,3),'k');
%axis([1 4000 0 7e-7])
xlabel('# of steps')
ylabel('backward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.4 unbounded','r = 0.1 spectral clustering')


%%
load('flux_avg_1ball_1step_beta_20_sc_200_2000_20.txt')
load('flux_avg_1ball_1step_beta_20_sc_200_4000_20.txt')
load('flux_avg_1ball_1step_beta_20_sc_400_8000_20.txt')
load('flux_avg_1ball_1step_beta_20_sc_600_12000_20.txt')
load('flux_avg_1ball_1step_beta_20_sc_600_12000_200.txt')
load('flux_avg_1ball_1step_beta_20_sc_600_12000_100_v2.txt')
load('flux_avg_1ball_1step_beta_20_sc_600_12000_200_v2.txt')
load('flux_avg_1ball_1step_beta_20_sc_600_12000_400_v2.txt')
flux_avg_1ball_1step_beta_20_r_0_1 = load('flux_avg_1ball_1step_beta_20_r_0.1.txt');
flux_avg_1ball_1step_beta_20_r_0_2 = load('flux_avg_1ball_1step_beta_20_r_0.2.txt');
flux_avg_1ball_1step_beta_20_r_0_4 = load('flux_avg_1ball_1step_beta_20_r_0.4.txt');
flux_avg_1ball_1step_beta_20_r_0_8 = load('flux_avg_1ball_1step_beta_20_r_0.8.txt');
flux_avg_1ball_1step_beta_20_r_0_1_v2 = load('flux_avg_1ball_1step_beta_20_r_0.1_v2.txt');
flux_avg_1ball_1step_beta_20_r_0_8_v2 = load('flux_avg_1ball_1step_beta_20_r_0.8_v2.txt');

load('flux_std_1ball_1step_beta_20_sc_200_2000_20.txt')
load('flux_std_1ball_1step_beta_20_sc_200_4000_20.txt')
load('flux_std_1ball_1step_beta_20_sc_400_8000_20.txt')
load('flux_std_1ball_1step_beta_20_sc_600_12000_20.txt')
load('flux_std_1ball_1step_beta_20_sc_600_12000_200.txt')
load('flux_std_1ball_1step_beta_20_sc_600_12000_100_v2.txt')
load('flux_std_1ball_1step_beta_20_sc_600_12000_200_v2.txt')
load('flux_std_1ball_1step_beta_20_sc_600_12000_400_v2.txt')
flux_std_1ball_1step_beta_20_r_0_1 = load('flux_std_1ball_1step_beta_20_r_0.1.txt');
flux_std_1ball_1step_beta_20_r_0_2 = load('flux_std_1ball_1step_beta_20_r_0.2.txt');
flux_std_1ball_1step_beta_20_r_0_4 = load('flux_std_1ball_1step_beta_20_r_0.4.txt');
flux_std_1ball_1step_beta_20_r_0_8 = load('flux_std_1ball_1step_beta_20_r_0.8.txt');
flux_std_1ball_1step_beta_20_r_0_1_v2 = load('flux_std_1ball_1step_beta_20_r_0.1_v2.txt');
flux_std_1ball_1step_beta_20_r_0_8_v2 = load('flux_std_1ball_1step_beta_20_r_0.8_v2.txt');


%%
figure; hold on;
freq = 200;
plot((1:20000),ones(20000,1)*4.8006567250e-08,'-.k');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_1(1:freq:660,1),flux_avg_1ball_1step_beta_20_r_0_1(1:freq:660,2),flux_std_1ball_1step_beta_20_r_0_1(1:freq:660,2),'r');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_1_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_1_v2(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_1_v2(1:freq:end,2),'r');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_2(1:freq:2229,1),flux_avg_1ball_1step_beta_20_r_0_2(1:freq:2229,2),flux_std_1ball_1step_beta_20_r_0_2(1:freq:2229,2),'g');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_4(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_8(1:freq:end,2),'m');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_8_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8_v2(1:freq:end,2),flux_std_1ball_1step_beta_20_r_0_8_v2(1:freq:end,2),'m');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,2),'r');
%errorbar(flux_avg_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,2),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,2),'r');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,2),'g');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,2),flux_std_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,2),'b');
%axis([1 7000 0 4e-6])
xlabel('# of steps')
ylabel('forward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.8 unbounded','r = 0.1 spectral clustering')


%%
figure; hold on;
freq = 200;
plot((1:20000),ones(20000,1)*4.8006567250e-08,'-.k');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_1(1:freq:660,1),flux_avg_1ball_1step_beta_20_r_0_1(1:freq:660,3),flux_std_1ball_1step_beta_20_r_0_1(1:freq:660,3),'r');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_1_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_1_v2(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_1_v2(1:freq:end,3),'r');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_2(1:freq:2229,1),flux_avg_1ball_1step_beta_20_r_0_2(1:freq:2229,3),flux_std_1ball_1step_beta_20_r_0_2(1:freq:2229,3),'g');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_4(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_4(1:freq:end,3),'b');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_8(1:freq:end,3),'m');
%errorbar(flux_avg_1ball_1step_beta_20_r_0_8_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_r_0_8_v2(1:freq:end,3),flux_std_1ball_1step_beta_20_r_0_8_v2(1:freq:end,3),'m');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_200_2000_20(1:freq:end,3),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_200_4000_20(1:freq:end,3),'r');
%errorbar(flux_avg_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_400_8000_20(1:freq:end,3),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_600_12000_20(1:freq:end,3),'b');
%errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_600_12000_200(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_600_12000_100_v2(1:freq:end,3),'r');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_600_12000_200_v2(1:freq:end,3),'g');
errorbar(flux_avg_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,1),flux_avg_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,3),flux_std_1ball_1step_beta_20_sc_600_12000_400_v2(1:freq:end,3),'b');
%axis([1 7000 0 10e-7])
xlabel('# of steps')
ylabel('backward flux')
%legend('exact','r = 0.1 unbounded','r = 0.2 unbounded','r = 0.8 unbounded','r = 0.1 spectral clustering')

