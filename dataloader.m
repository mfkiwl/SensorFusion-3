function [bag_data,rlog_data] = dataloader(date)
% dataloader(date) loads bag_data.mat, rlog_data.mat from input date
% Implemented by JinHwan Jeon, 2023
    rlog_data_path = strcat('/media/jinhwan/JinHwan/SJ_Dataset/2023/',date,'/rlog/rlog_data.mat');
    bag_data_path = strcat('/media/jinhwan/JinHwan/SJ_Dataset/2023/',date,'/bag/bag_data.mat');
    rlog_data = load(rlog_data_path);
    bag_data = load(bag_data_path,'output');
end