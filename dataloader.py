from math import inf
import pickle
import pandas as pd
import pymap3d as pm
import numpy as np
from numpy.linalg import norm, inv
from scipy import interpolate
from bisect import bisect_left
import matplotlib.pyplot as plt
import scipy.io as sio
import os

class DataLoader:
    """
    Load and extract data into usable format from .pkl files
    Need to change property 'names' to change files to read

    Original version interpolates timestamps but found some error!
    Code is modified to perform simple data extraction only!

    Implemented by JinHwan Jeon, 2022
    """
    def __init__(self, file_path, save_path):
        self.FOI = file_path
        self.save_path = save_path

        # Files to read: add more if needed
        self.names = ['carState.pkl',
                      'gpsLocationExternal.pkl',
                      'liveLocationKalman.pkl',
                      'modelV2.pkl',
                      'sensorEvents.pkl',
                      'snapRaw.pkl']

        # snapRaw.pkl
        # snapMatched.pkl

        self.raw_data = {} # Dictionary to save loaded data files (Raw)
        self.data = {} # Semi-processed Data
        self.proc_data = {} # Processed Data
        self.unixtime_base = 0 # Local time base

    def load(self):
        print("[Reading pkl files...]")
        for name in self.names:
            with open(self.FOI + name, 'rb') as f:
                self.raw_data[name[0:-4]] = pickle.load(f)
        
        # Extract IMU readings 
        self.extractData()

        self.save2mat()

    def extractData(self):
        
        # Read CAN data
        DOI = self.raw_data['carState']
        
        left = DOI['leftBlinker'].to_numpy()
        right = DOI['rightBlinker'].to_numpy()
        whl_spd = DOI['wheelSpeeds.rl'].to_numpy()
        
        can_t = DOI.index

        self.data['can'] = {}
        self.data['can']['whl_spd'] = []
        self.data['can']['t'] = []
        self.data['can']['leftBlinker'] = []
        self.data['can']['rightBlinker'] = []

        for i in range(len(can_t)):
            self.data['can']['whl_spd'].append(whl_spd[i])
            self.data['can']['t'].append(can_t[i])
            self.data['can']['leftBlinker'].append(left[i])
            self.data['can']['rightBlinker'].append(right[i])

        # Read GNSS data

        DOI = self.raw_data['gpsLocationExternal']
        DOI['idx'] = list(range(len(DOI.index)))

        # DOI = DOI.loc[DOI['accuracy'] < 5] # Discard data with horizontal accuracy higher than 5m 

        lat = DOI['latitude'].to_numpy()
        lon = DOI['longitude'].to_numpy()
        alt = DOI['altitude'].to_numpy()
        hAcc = DOI['accuracy'].to_numpy()
        vAcc = DOI['verticalAccuracy'].to_numpy()
        speed = DOI['speed'].to_numpy()
        vNED = DOI['vNED'].to_numpy()
        bearing = DOI['bearingDeg'].to_numpy()
        bearingAcc = DOI['bearingAccuracyDeg'].to_numpy()
        raw_idx = np.array(DOI['idx'])

        self.data['gnss'] = {}
        self.data['gnss']['t'] = []
        self.data['gnss']['t_bias'] = [] # 
        self.data['gnss']['pos'] = []
        self.data['gnss']['lla_ref'] = [36.5222,  127.3032,  181.1957] # Reference for PPK
        self.data['gnss']['hAcc'] = []
        self.data['gnss']['vAcc'] = []
        self.data['gnss']['vNED'] = []
        self.data['gnss']['bearing'] = []
        self.data['gnss']['bearingAcc'] = []
        self.data['gnss']['raw_idx'] = []
        self.data['gnss']['speed'] = []

        intv_cnt = 0
        static_time_bias = [0] # time value to be subtracted for data continuity
        gnss_t = DOI.index

        # Clustering GNSS data for non-static intervals
        for i in range(len(gnss_t)):
            self.data['gnss']['t'].append(gnss_t[i])
            self.data['gnss']['pos'].append([lat[i],lon[i],alt[i]])
            self.data['gnss']['hAcc'].append(hAcc[i])
            self.data['gnss']['vAcc'].append(vAcc[i])
            self.data['gnss']['vNED'].append([val for val in vNED[i]])
            self.data['gnss']['bearing'].append(bearing[i])
            self.data['gnss']['bearingAcc'].append(bearingAcc[i])
            self.data['gnss']['raw_idx'].append(raw_idx[i])
            self.data['gnss']['speed'].append(speed[i])

        
        lat = [x[0] for x in self.data['gnss']['pos']]
        lon = [x[1] for x in self.data['gnss']['pos']]

        plt.plot(lon,lat,'r.')
        plt.show()
        plt.axis('equal')

        # Read and classify IMU data
        
        DOI = self.raw_data['sensorEvents']['']
        self.data['imu'] = {}

        IMU = {}
        self.data['imu']['accel'] = {}
        self.data['imu']['accel']['t'] = []
        self.data['imu']['accel']['meas'] = []
        self.data['imu']['gyro'] = {}
        self.data['imu']['gyro']['t'] = []
        self.data['imu']['gyro']['meas'] = []

        cnt = 0
        for row in DOI:
            for el in row:
                if 'gyroUncalibrated' in el.keys():
                    self.data['imu']['gyro']['t'].append(DOI.index[cnt])
                    meas = el['gyroUncalibrated']['v']
                    self.data['imu']['gyro']['meas'].append([meas[0],meas[1],meas[2]])

                elif 'acceleration' in el.keys():
                    self.data['imu']['accel']['t'].append(DOI.index[cnt])
                    meas = el['acceleration']['v']
                    self.data['imu']['accel']['meas'].append([meas[0],meas[1],meas[2]])

            cnt += 1


        # ax = []; ay = []; az = []; wx = []; wy = []; wz = []; t = []
        # for acc in self.data['imu']['accel']['meas']:
        #     ax.append(acc[0])
        #     ay.append(acc[1])
        #     az.append(acc[2])
        # for gyro in self.data['imu']['gyro']['meas']:
        #     wx.append(gyro[0])
        #     wy.append(gyro[1])
        #     wz.append(gyro[2])

        # plt.figure(1)
        # plt.plot(ax,'r')
        # plt.plot(ay,'g')
        # plt.plot(az,'b')
        # plt.figure(2)
        # plt.plot(wx,'r')
        # plt.plot(wy,'g')
        # plt.plot(wz,'b')

        # plt.show()

        """
        * Read lane measurements
        No need to classify (matching will be done)
        Just save raw data into interpretable form

        """
        DOI = self.raw_data['modelV2']
        self.data['lane'] = {}
        self.data['lane']['t'] = []
        self.data['lane']['x'] = DOI['laneLines.0.x'].to_numpy()[0]
        self.data['lane']['x_inter'] = np.linspace(0,100,11)

        lly = DOI['laneLines.0.y'].to_numpy(); llystd = DOI['laneLines.0.yStd'].to_numpy()
        ly = DOI['laneLines.1.y'].to_numpy(); lystd = DOI['laneLines.1.yStd'].to_numpy()
        ry = DOI['laneLines.2.y'].to_numpy(); rystd = DOI['laneLines.2.yStd'].to_numpy()
        rry = DOI['laneLines.3.y'].to_numpy(); rrystd = DOI['laneLines.3.yStd'].to_numpy()
        ley = DOI['roadEdges.0.y'].to_numpy(); leystd = DOI['roadEdges.0.yStd'].to_numpy()
        rey = DOI['roadEdges.1.y'].to_numpy(); reystd = DOI['roadEdges.1.yStd'].to_numpy() 

        llz = DOI['laneLines.0.z'].to_numpy(); llzstd = DOI['laneLines.0.zStd'].to_numpy()
        lz = DOI['laneLines.1.z'].to_numpy(); lzstd = DOI['laneLines.1.zStd'].to_numpy()
        rz = DOI['laneLines.2.z'].to_numpy(); rzstd = DOI['laneLines.2.zStd'].to_numpy()
        rrz = DOI['laneLines.3.z'].to_numpy(); rrzstd = DOI['laneLines.3.zStd'].to_numpy()
        lez = DOI['roadEdges.0.z'].to_numpy(); lezstd = DOI['roadEdges.0.zStd'].to_numpy()
        rez = DOI['roadEdges.1.z'].to_numpy(); rezstd = DOI['roadEdges.1.zStd'].to_numpy()

        prob = DOI['laneLineProbs'].to_numpy()
        
        self.data['lane']['lly'] = []; self.data['lane']['llystd'] = []
        self.data['lane']['ly'] = []; self.data['lane']['lystd'] = []
        self.data['lane']['ry'] = []; self.data['lane']['rystd'] = []
        self.data['lane']['rry'] = []; self.data['lane']['rrystd'] = []
        self.data['lane']['ley'] = []; self.data['lane']['leystd'] = []
        self.data['lane']['rey'] = []; self.data['lane']['reystd'] = []

        self.data['lane']['llz'] = []; self.data['lane']['llzstd'] = []
        self.data['lane']['lz'] = []; self.data['lane']['lzstd'] = []
        self.data['lane']['rz'] = []; self.data['lane']['rzstd'] = []
        self.data['lane']['rrz'] = []; self.data['lane']['rrzstd'] = []
        self.data['lane']['lez'] = []; self.data['lane']['lezstd'] = []
        self.data['lane']['rez'] = []; self.data['lane']['rezstd'] = []

        self.data['lane']['prob'] = []

        for i in range(len(DOI.index)):

            self.data['lane']['t'].append(DOI.index[i])
            
            # 1D Interpolation to prevent overfitting at near preview distances
            # Y
            self.data['lane']['lly'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in lly[i]]))
            self.data['lane']['llystd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in llystd[i]]))
            self.data['lane']['ly'].append(np.interp(self.data['lane']['x_inter'],
                                                     self.data['lane']['x'],
                                                     [-val for val in ly[i]]))
            self.data['lane']['lystd'].append(np.interp(self.data['lane']['x_inter'],
                                                        self.data['lane']['x'],
                                                        [val for val in lystd[i]]))
            self.data['lane']['ry'].append(np.interp(self.data['lane']['x_inter'],
                                                     self.data['lane']['x'],
                                                     [-val for val in ry[i]]))
            self.data['lane']['rystd'].append(np.interp(self.data['lane']['x_inter'],
                                                        self.data['lane']['x'],
                                                        [val for val in rystd[i]]))
            self.data['lane']['rry'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in rry[i]]))
            self.data['lane']['rrystd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in rrystd[i]]))

            self.data['lane']['ley'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in ley[i]]))
            self.data['lane']['leystd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in leystd[i]]))
            self.data['lane']['rey'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in rey[i]]))
            self.data['lane']['reystd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in reystd[i]]))

            # Z
            self.data['lane']['llz'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in llz[i]]))
            self.data['lane']['llzstd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in llzstd[i]]))
            self.data['lane']['lz'].append(np.interp(self.data['lane']['x_inter'],
                                                     self.data['lane']['x'],
                                                     [-val for val in lz[i]]))
            self.data['lane']['lzstd'].append(np.interp(self.data['lane']['x_inter'],
                                                        self.data['lane']['x'],
                                                        [val for val in lzstd[i]]))
            self.data['lane']['rz'].append(np.interp(self.data['lane']['x_inter'],
                                                     self.data['lane']['x'],
                                                     [-val for val in rz[i]]))
            self.data['lane']['rzstd'].append(np.interp(self.data['lane']['x_inter'],
                                                        self.data['lane']['x'],
                                                        [val for val in rzstd[i]]))
            self.data['lane']['rrz'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in rrz[i]]))
            self.data['lane']['rrzstd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in rrzstd[i]]))

            self.data['lane']['lez'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in lez[i]]))
            self.data['lane']['lezstd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in lezstd[i]]))
            self.data['lane']['rez'].append(np.interp(self.data['lane']['x_inter'],
                                                      self.data['lane']['x'],
                                                      [-val for val in rez[i]]))
            self.data['lane']['rezstd'].append(np.interp(self.data['lane']['x_inter'],
                                                         self.data['lane']['x'],
                                                         [val for val in rezstd[i]]))

            self.data['lane']['prob'].append([val for val in prob[i]])

        DOI = self.raw_data['snapRaw']
        t = DOI['timestamp_snap'].to_numpy()
        lat = DOI['latitude_snap'].to_numpy()
        lon = DOI['longitude_snap'].to_numpy()
        dist = DOI['distance_snap'].to_numpy()

        self.data['snap_raw'] = {}
        self.data['snap_raw']['t'] = t
        self.data['snap_raw']['lat'] = lat
        self.data['snap_raw']['lon'] = lon
        self.data['snap_raw']['dist'] = dist

        # DOI = self.raw_data['snapMatched']
        # t = DOI['timestamp_snap'].to_numpy()
        # lat = DOI['latitude_snap'].to_numpy()
        # lon = DOI['longitude_snap'].to_numpy()
        # dist = DOI['distance_snap'].to_numpy()

        # self.data['snap'] = {}
        # self.data['snap']['t'] = t
        # self.data['snap']['lat'] = lat
        # self.data['snap']['lon'] = lon
        # self.data['snap']['dist'] = dist

    def save2mat(self):
        """
        Save Dictonaries into .mat files
        
        """
        print("[Saving Files to MATLAB format...]")

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)
            
        sio.savemat(self.save_path+'gnss.mat',self.data['gnss'])
        sio.savemat(self.save_path+'lane.mat',self.data['lane'])
        sio.savemat(self.save_path+'can.mat',self.data['can'])
        sio.savemat(self.save_path+'imu.mat',self.data['imu'])
        sio.savemat(self.save_path+'snap_raw.mat',self.data['snap_raw'])
        # sio.savemat(self.save_path+'snap.mat',self.data['snap'])
        
def findIntvs(idxs,t):
    n = len(idxs)
    intvs = []
    lb = idxs[0]

    for i in range(n-1):
        if idxs[i+1] - idxs[i] != 1:
            ub = idxs[i]
            intvs.append([t[lb-1],t[ub+1]])
                
            lb = idxs[i+1]
        if i == n-2:
            ub = idxs[-1]
            intvs.append([t[lb-1], t[ub+1]])

    return intvs

def closestIdx(data,val):
    """
    Assumes data is sorted. Returns closest data index.
    If two numbers are equally close, return the smallest number.

    Algorithm complexity is O(log n), much better than simply using "argmin" for large dataset
    """
    pos = bisect_left(data,val)
    if pos == 0:
        return 0
    if pos == len(data):
        return len(data)-1
    before = data[pos - 1]
    after = data[pos]
    if after - val < val - before:
        return pos
    else:
        return pos-1
