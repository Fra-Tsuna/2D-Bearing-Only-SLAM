close all;
clear;
clc;

addpath "./g2o_wrapper";
addpath "./dataset";
source "functions.m";

##Extraction of data from .g2o files

[land_gt, pose_gt, trans_gt, obs_gt] = loadG2o('slam2D_bearing_only_ground_truth.g2o');
[land_ig, pose_ig, trans_ig, obs_ig] = loadG2o('slam2D_bearing_only_initial_guess.g2o');

land_gt=land_gt(2:end);
pose_gt=pose_gt(2:end);
trans_gt=trans_gt(2:end);
obs_gt=obs_gt(2:end);

land_ig=land_ig(2:end);
pose_ig=pose_ig(2:end);
trans_ig=trans_ig(2:end);
obs_ig=obs_ig(2:end);


#init poses (gt and ig)
Nr=length(pose_gt);

Xr_gt = zeros(3,3,Nr);
Xr_ig = zeros(3,3,Nr);

pose_id2index_gt = zeros(1,Nr); 
pose_id2index_ig = zeros(1,Nr); 

for i=1:Nr
    x_gt=[pose_gt(i).x, pose_gt(i).y, pose_gt(i).theta];
    Xr_gt(:,:,i) = v2t(x_gt);
    pose_id2index_gt(i)= pose_gt(i).id;

    x_ig=[pose_ig(i).x, pose_ig(i).y, pose_ig(i).theta];
    Xr_ig(:,:,i) = v2t(x_ig);
    pose_id2index_ig(i)= pose_ig(i).id;
end

#init landmarks positions (gt)
Nl=length(land_gt);
Xl = zeros(2,Nl);
land_id2index = zeros(1,Nl); 

for i=1:Nl
    Xl(:,i) = [land_gt(i).x_pose; land_gt(i).y_pose];
    land_id2index(i) = land_gt(i).id;
end

state_dim = 3*Nr + 2*Nl;

#find total number of measurements
Ml_ig=0;
for i=1:length(obs_ig)
    m=length(obs_ig(i).observation);
    Ml_ig += m;
end

#init bearing measurements 
Zl_ig=zeros(1,Ml_ig);
associations_Zl_ig = zeros(2,Ml_ig);

#create for each bearing measurement an association between the pose that recorded it
#and the landmark that originated it while extracting them
m=1;
for i=1:length(obs_ig)
    observations = obs_ig(i);
    pose_id = observations.pose_id;
    pose_index = find(pose_id2index_ig==pose_id);
    num_observations = length(observations.observation);
    for j=1:num_observations
        measurement = observations.observation(j);
        land_id = measurement.id;
        land_index = find(land_id2index==land_id);
        Zl_ig(m)=measurement.bearing;
        associations_Zl_ig(:,m) = [pose_index; land_index];
        m++;
    end
end

#init odometry measurements
Mr_ig = length(trans_ig);
Zr_ig = zeros(3,3,Mr_ig);
associations_Zr_ig = zeros(2,Mr_ig);


#create for each odometry measurement an association between the previous pose and
#the subsequent pose
m=1;
for i=1:Mr_ig
    transition = trans_ig(i);
    pose_from_id = transition.id_from;
    pose_to_id = transition.id_to;
    pose_from_index = find(pose_id2index_ig==pose_from_id);
    pose_to_index = find(pose_id2index_ig==pose_to_id);
    associations_Zr_ig(:,i) = [pose_from_index; pose_to_index];
    Zr_ig(:,:,i) = v2t(transition.v);
end

#landmarks triangulation for landmarks initial guess position


figure(1);
grid;
plot(Xl(1,:),Xl(2,:),'b*',"linewidth",2);
hold on;
plot(squeeze(Xr_gt(1,3,:)),squeeze(Xr_gt(2,3,:)),'b*-',"linewidth",2);
hold on;
plot(squeeze(Xr_ig(1,3,:)),squeeze(Xr_ig(2,3,:)),'g*-',"linewidth",2);