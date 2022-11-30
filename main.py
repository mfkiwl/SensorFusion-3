from dataloader import DataLoader
import matplotlib.pyplot as plt
import numpy as np

"""
Main File for data processing and optimization
Implemented by JinHwan Jeon, 2022
"""
# Load Data (.pkl files)
# path = 'D:/SJ_Dataset/2022-08-05/2022-08-05--02-28-29/'
# path = 'D:/SJ_Dataset/20220805_022829/'
# file_path = 'C:/Users/user/Desktop/2022/연구/SejongDataset/2022-08-05--02-41-03/'
# save_path = 'C:/Users/user/Desktop/2022/연구/SejongDataset/2022-08-05--02-41-03/results/'

# file_path = 'C:/Users/user/Desktop/2022/연구/SejongDataset/2022-08-05--02-57-24/'
# save_path = 'C:/Users/user/Desktop/2022/연구/SejongDataset/2022-08-05--02-57-24/results/'

# 'D:/SJ_Dataset/2022-08-05/2022-08-05--02-41-03/'
# 'D:/SJ_Dataset/2022-08-05/2022-08-05--02-57-24/'
# 'D:/SJ_Dataset/2022-08-05/2022-08-05--03-13-16/'
# 'D:/SJ_Dataset/2022-08-05/2022-08-05--03-45-16/'
# 'D:/SJ_Dataset/2022-08-05/2022-08-05--04-19-33/'
file_path = 'D:/SJ_Dataset/2022-08-05/2022-08-05--02-41-03/'
save_path = 'D:/SJ_Dataset/2022-08-05/2022-08-05--02-41-03/results/'

# path = ''
dataset = DataLoader(file_path,save_path)
dataset.load()

"""
Check acceleration values to match axis correctly

For other dataset with multiple static scenes, check if time_bias terms are removed correctly

"""
# ax = []; ay = []; wz = []; t = []
# for imuReadings in dataset.proc_data['imu']:
#     # print(imuReadings)
#     acc = imuReadings['accel']
#     for row in acc:
#         ax.append(row[0])
#         ay.append(row[1])
    
#     gyro = imuReadings['gyro']
#     for row in gyro:
#         wz.append(row[2])

#     t.extend([val for val in imuReadings['t']])

# plt.figure(4)
# plt.plot(ax,'r')
# plt.plot(ay,'b')
# plt.plot(wz,'g')

# plt.figure(5)
# plt.plot(t,'k.')
# plt.show()



