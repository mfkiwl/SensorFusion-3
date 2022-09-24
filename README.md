# SensorFusion
SensorFusion using MATLAB

## Target
Vehicle attitude, trajectory and lane map reconstruction for sparse feature, GNSS degraded, high speed drive environment 

## Files
### Starting with raw data

- To start with raw data, run 'main.py' file first with appropriate file path (raw data files are not included in this repo)

- 'dataloader.py' file reads '.pkl' format files and extracts data into dictionary variable, and finally converts to '.mat' file

- Modify 'dataloader.py' to extract different raw data from other '.pkl' files

### Starting with extracted data
#### 'dataprocessor.m'

- After running 'dataloader.py', there will be various '.mat' format data files

- Running one section in 'main.m' will read the data and process through them for creating dataset of directly usable format

- Timestamps are interpolated, GNSS measurements with poor accuracy are filtered for optimization stability

- Finally, IMU readings are clustered for easier preintegration

#### 'optimizer.m'

- Running the remaining sections of 'main.m' will use 'optimizer.m', which is the main part of this research

- Using Sparse Non-Linear Least Squares optimization algorithms, 'optimizer.m' solves Sensor Fusion problem

- There are currently 3 modes possible, 'basic', 'partial', '2-phase'. These modes are classified based on the sensor data used for optimization

- Read 'optimizer.m' for more detailed explanation