function output = rosbag_convert(filename)
% ROSBAG_CONVERT(filename) is a (.bag) file to (.mat) data converter.
% Edit or add topic name to read messages from other topics

    bag = rosbag(filename);
    output = struct();
    
    %% 
    
    temptopic = "/RT/gps/vel";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    templen = length(tempmsg);
    msgname = 'rtvel';
    output.(msgname) = struct();
    output.(msgname).t = zeros(1,templen);
    output.(msgname).vx = zeros(1,templen);
    output.(msgname).vy = zeros(1,templen);
    output.(msgname).vz = zeros(1,templen);
    
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
        output.(msgname).vx(i) = msg.Twist.Twist.Linear.X;
        output.(msgname).vy(i) = msg.Twist.Twist.Linear.Y;
        output.(msgname).vz(i) = msg.Twist.Twist.Linear.Z;
    end
    
    %% 
    
    temptopic = "/RT/gps/fix";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    templen = length(tempmsg);
    msgname = 'rtgps';
    output.(msgname) = struct();
    output.(msgname).t = zeros(1,templen);
    output.(msgname).lat = zeros(1,templen);
    output.(msgname).lon = zeros(1,templen);
    output.(msgname).alt = zeros(1,templen);
    
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
        output.(msgname).lat(i) = msg.Latitude;
        output.(msgname).lon(i) = msg.Longitude;
        output.(msgname).alt(i) = msg.Altitude;
    end
    
    %%
    
    temptopic = "/RT/imu/data";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    templen = length(tempmsg);
    msgname = 'rtimu';
    output.(msgname) = struct();
    output.(msgname).t = zeros(1,templen);
    output.(msgname).ax = zeros(1,templen);
    output.(msgname).ay = zeros(1,templen);
    output.(msgname).az = zeros(1,templen);
    output.(msgname).wx = zeros(1,templen);
    output.(msgname).wy = zeros(1,templen);
    output.(msgname).wz = zeros(1,templen);
    output.(msgname).rx = zeros(1,templen);
    output.(msgname).ry = zeros(1,templen);
    output.(msgname).rz = zeros(1,templen);
    
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
        output.(msgname).ax(i) = msg.LinearAcceleration.X;
        output.(msgname).ay(i) = -msg.LinearAcceleration.Y;
        output.(msgname).az(i) = msg.LinearAcceleration.Z;
        output.(msgname).wx(i) = msg.AngularVelocity.X;
        output.(msgname).wy(i) = -msg.AngularVelocity.Y;
        output.(msgname).wz(i) = msg.AngularVelocity.Z;
        temp = quat2eul([msg.Orientation.W, msg.Orientation.X, msg.Orientation.Y, msg.Orientation.Z]);
        output.(msgname).rx(i) = temp(3);
        output.(msgname).ry(i) = -temp(2);
        output.(msgname).rz(i) = temp(1);
    end
    
    % %%
    % 
    % temptopic = "/openpilot/carState"
    % tempbag = select(bag, 'Topic' , temptopic);
    % tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    % templen = length(tempmsg);
    % msgname = 'opcarstate';
    % output.(msgname) = struct();
    % output.(msgname).t = zeros(1,templen);
    % output.(msgname).delta = zeros(1,templen);
    % output.(msgname).wsf = zeros(1,templen);
    % output.(msgname).wsr = zeros(1,templen);
    % output.(msgname).v = zeros(1,templen);
    % for i = 1:length(tempmsg)
    %     msg = tempmsg{i};
    %     output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
    %     output.(msgname).delta(i) = msg.SteeringAngle.Data;
    %     output.(msgname).wsf(i) = (msg.WheelSpeeds.Fl.Data+msg.WheelSpeeds.Fr.Data)/2;
    %     output.(msgname).wsr(i) = (msg.WheelSpeeds.Rl.Data+msg.WheelSpeeds.Rr.Data)/2;
    %     output.(msgname).v(i) = msg.VEgo.Data;
    % end
    % 
    %%
    
    temptopic = "/filter/quaternion";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    templen = length(tempmsg);
    msgname = 'xs';
    output.(msgname) = struct();
    output.(msgname).t = zeros(1,templen);
    output.(msgname).rx = zeros(1,templen);
    output.(msgname).ry = zeros(1,templen);
    output.(msgname).rz = zeros(1,templen);
    output.(msgname).vx = zeros(1,templen);
    output.(msgname).vy = zeros(1,templen);
    output.(msgname).vz = zeros(1,templen);
    output.(msgname).ax = zeros(1,templen);
    output.(msgname).ay = zeros(1,templen);
    output.(msgname).az = zeros(1,templen);
    output.(msgname).wx = zeros(1,templen);
    output.(msgname).wy = zeros(1,templen);
    output.(msgname).wz = zeros(1,templen);
    
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
        temp = quat2eul([msg.Quaternion.W, msg.Quaternion.X, msg.Quaternion.Y, msg.Quaternion.Z]);
        output.(msgname).rx(i) = temp(3);
        output.(msgname).ry(i) = temp(2);
        output.(msgname).rz(i) = temp(1);
    end
    
    temptopic = "/filter/velocity";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).vx(i) = msg.Vector.X;
        output.(msgname).vy(i) = msg.Vector.Y;
        output.(msgname).vz(i) = msg.Vector.Z;
    end
    
    temptopic = "/imu/acceleration";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).ax(i) = msg.Vector.X;
        output.(msgname).ay(i) = msg.Vector.Y;
        output.(msgname).az(i) = msg.Vector.Z;
    end
    
    temptopic = "/imu/angular_velocity";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).wx(i) = msg.Vector.X;
        output.(msgname).wy(i) = msg.Vector.Y;
        output.(msgname).wz(i) = msg.Vector.Z;
    end
    
    temptopic = "/imu/data";
    tempbag = select(bag, 'Topic' , temptopic);
    tempmsg = readMessages(tempbag, 'DataFormat', 'struct');
    
    output.(msgname).wcov = zeros(9,templen);
    output.(msgname).acov = zeros(9,templen);
    
    for i = 1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).wcov(:,i) = msg.AngularVelocityCovariance;
        output.(msgname).acov(:,i) = msg.LinearAccelerationCovariance;
    end
    
    %% 
    temptopic = "/gnss";
    tempbag = select(bag, 'Topic', temptopic);
    tempmsg = readMessages(tempbag,'DataFormat','struct');
    templen = length(tempmsg);
    
    msgname = 'rtk';
    output.(msgname) = struct();
    output.(msgname).t = zeros(1,templen);
    output.(msgname).lat = zeros(1,templen);
    output.(msgname).lon = zeros(1,templen);
    output.(msgname).alt = zeros(1,templen);
    output.(msgname).cov = zeros(9,templen);

    for i=1:length(tempmsg)
        msg = tempmsg{i};
        output.(msgname).t(i) = double(msg.Header.Stamp.Sec) + double(msg.Header.Stamp.Nsec)*1e-9;
        output.(msgname).lat(i) = msg.Latitude;
        output.(msgname).lon(i) = msg.Longitude;
        output.(msgname).alt(i) = msg.Altitude;
        output.(msgname).cov(:,i) = msg.PositionCovariance;
    end

end