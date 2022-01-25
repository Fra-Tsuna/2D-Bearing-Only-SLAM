close all;
clear;
clc;

addpath "./g2o_wrapper";
addpath "./dataset";
source "functions.m";
source "do_Least_Square_Magic.m"

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
Xl_gt = zeros(2,Nl);
land_id2index = zeros(1,Nl); 

for i=1:Nl
    Xl_gt(:,i) = [land_gt(i).x_pose; land_gt(i).y_pose];
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
Xl_ig =zeros(2,Nl);

[associations_sorted,indices] = sort(associations_Zl_ig(2,:));
landmark_ids = unique(associations_sorted);
obs_per_landmarks = zeros(size(landmark_ids));

for i=1:length(landmark_ids)
    landmark_id = landmark_ids(i);
    bearing_values = Zl_ig(find(associations_Zl_ig(2,:)==landmark_id));
    rob_poses_from_land_id = associations_Zl_ig(1,:)(find(associations_Zl_ig(2,:)==landmark_id));
    rob_poses = zeros(3,length(bearing_values));
    for j = 1:length(bearing_values)
        rob_poses(:,j) = t2v(Xr_ig(:,:,rob_poses_from_land_id(j)));
        bearing_values(j) = bearing_values(j) + rob_poses(3,j);
        bearing_values(j) = atan2(sin(bearing_values(j)),cos(bearing_values(j)));
    end

    obs_per_landmarks(landmark_id) = sum(associations_sorted==landmark_id);

    %%now we look for the 2 best poses in order to do the triangulation
    A = zeros(obs_per_landmarks(landmark_id),obs_per_landmarks(landmark_id));
    for m=1:obs_per_landmarks(landmark_id)-1
        for n=(m+1):obs_per_landmarks(landmark_id)
            A(m,n) = box_minus(bearing_values(m),bearing_values(n));
            if A(m,n)>pi/2
                A(m,n)=pi-A(m,n);
            endif
        end
    end
    [rob_pose_1, rob_pose_2] = find(A==max(max(A)));

    if(obs_per_landmarks(landmark_id)==1)
        Xl_ig(:,landmark_id) = rob_poses(1:2,1)+[cos(rob_poses(3,1)); sin(rob_poses(3,1))];
    else
        p11 = rob_poses(1:2,rob_pose_1);
        theta1 = bearing_values(rob_pose_1);
        p12 = p11 + [cos(theta1); sin(theta1)];
        p21 = rob_poses(1:2,rob_pose_2);
        theta2 = bearing_values(rob_pose_2);
        p22 = p21 + [cos(theta2); sin(theta2)];
        Xl_ig(:,landmark_id) = intersect_lines(p11,p12,p21,p22);
    endif
end

iterations=100;

[Xr_corr, Xl_corr, chi_stats_poses, chi_stats_landmarks] = ...
    do_Least_Square_Magic(Xr_ig, Xl_ig, Zl_ig, associations_Zl_ig, Zr_ig, associations_Zr_ig, iterations);


#plots and figures
figure(1);
plot(Xl_gt(1,:),Xl_gt(2,:),'b*',"linewidth",2);
hold on;
plot(squeeze(Xr_gt(1,3,:)),squeeze(Xr_gt(2,3,:)),'b*-',"linewidth",2);
hold on;
plot(Xl_ig(1,:),Xl_ig(2,:),'g*',"linewidth",2);
hold on;
plot(squeeze(Xr_ig(1,3,:)),squeeze(Xr_ig(2,3,:)),'g*-',"linewidth",2);
grid;