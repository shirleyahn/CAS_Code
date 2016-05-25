import numpy as np
import os

os.system('awk \'{print $5\'} time_record.txt > time1.txt')
os.system('awk \'{print $8\'} time_record.txt > time2.txt')
os.system('awk \'{print $11\'} time_record.txt > time3.txt')

time1 = np.loadtxt('time1.txt')
time2 = np.loadtxt('time2.txt')
time3 = np.loadtxt('time3.txt')

total_time = np.zeros((time1.shape[0],))
accumulation = 0.0
for i in range(time1.shape[0]):
    accumulation += time1[i]+time2[i]+time3[i]
    total_time[i] += accumulation

os.system('rm time1.txt')
os.system('rm time2.txt')
os.system('rm time3.txt')
os.system('rm time_record.txt')

np.savetxt('total_time.txt', total_time, fmt=' %+1.5f')