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
pose_id2index = zeros(Nr); 

for i=1:Nr
    x_gt=[pose_gt(i).x, pose_gt(i).y, pose_gt(i).theta];
    Xr_gt(:,:,i) = v2t(x_gt);
    pose_id2index(i)= pose_gt(i).id;
end

#init landmarks positions
Nl=length(land_gt);
Xl_gt = zeros(2,Nl);
land_id2index = zeros(Nl); 

for i=1:Nl
    Xl_gt(:,i) = [land_gt(i).x_pose; land_gt(i).y_pose];
    land_id2index(i) = land_gt(i).id;
end

state_dim = 3*Nr + 2*Nl;

#init measurements
Ml=0;
for i=1:length(obs_gt)
    m=length(obs_gt(i).observation);
    Ml += m;
end

Zl=zeros(Ml);
