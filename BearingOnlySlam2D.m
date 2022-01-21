close all;
clear;
clc;

addpath "./g2o_wrapper";
addpath "./dataset";
source "functions.m";

##Extraction of data from .g2o files
##Groundtruth

[land_gt, pose_gt, trans_gt, obs_gt] = loadG2o('slam2D_bearing_only_ground_truth.g2o');

#init poses
Nr=length(pose_gt);
Xr_gt = zeros(3,3,Nr);
pose_id2index = zeros(1,Nr); 

for i=1:Nr
    x_gt=[pose_gt(i).x, pose_gt(i).y, pose_gt(i).theta];
    Xr_gt(:,:,i) = v2t(x_gt);
    pose_id2index(i)= pose_gt(i).id;
end

#init landmarks positions
Nl=length(land_gt);
Xl_gt = zeros(2,Nl);
land_id2index = zeros(1,Nl); 

for i=1:Nl
    Xl_gt(:,i) = [land_gt(i).x_pose; land_gt(i).y_pose];
    land_id2index(i) = land_gt(i).id;
end

state_dim = 3*Nr + 2*Nl;

#find total number of measurements
Ml=0;
for i=1:length(obs_gt)
    m=length(obs_gt(i).observation);
    Ml += m;
end

#init bearing measurements
Zl=zeros(1,Ml);
associations_Zl = zeros(2,Ml);

#create for each bearing measurement an association between the pose that recorded it
#and the landmark that originated it while extracting them
m=1;
for i=1:length(obs_gt)
    observations = obs_gt(i);
    pose_id = observations.pose_id;
    pose_index = find(pose_id2index==pose_id);
    num_observations = length(observations.observation);
    for j=1:num_observations
        measurement = observations.observation(j);
        land_id = measurement.id;
        land_index = find(land_id2index==land_id);
        Zl(m)=measurement.bearing;
        associations_Zl(:,m) = [pose_index; land_index];
        m++;
    end
end

#init odometry measurements
Mr = length(trans_gt);
Zr = zeros(3,3,Mr);
associations_Zr = zeros(2,Mr);


#create for each odometry measurement an association between the previous pose and
#the subsequent pose
m=1;
for i=1:Mr
    transition = trans_gt(i);
    pose_from_id = transition.id_from;
    pose_to_id = transition.id_to;
    pose_from_index = find(pose_id2index==pose_from_id);
    pose_to_index = find(pose_id2index==pose_to_id);
    associations_Zr(:,i) = [pose_from_index; pose_to_index];
    Zr(:,:,i) = v2t(transition.v);
end