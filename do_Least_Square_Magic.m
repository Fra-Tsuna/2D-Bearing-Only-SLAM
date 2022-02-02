1;

function [Xr_corr, Xl_corr, chi_stats_poses, chi_stats_bearings] = do_Least_Square_Magic(Xr_ig, Xl_ig, Zl_ig, associations_Zl_ig, Zr_ig, associations_Zr_ig, iterations)
    Nr = size(Xr_ig)(3);
    Nl = size(Xl_ig)(2);
    Mr = size(Zr_ig)(3);
    Ml = size(Zl_ig)(2);


    state_dim = 3*Nr + 2*Nl;
    
    chi_stats_poses = zeros(1,iterations);
    chi_stats_bearings = zeros(1, iterations);

    Xr_corr = Xr_ig;
    Xl_corr = Xl_ig;

    disp("Starting Gauss-Newton algorithm...");

    for i=1:iterations
        printf("Iteration: %d\n", i);

        Hr = zeros(state_dim, state_dim);
        br = zeros(state_dim,1);
        %handle poses
        for m=1:Mr
            pose_i = associations_Zr_ig(1,m);
            pose_j = associations_Zr_ig(2,m);
            z = Zr_ig(:,:,m);
            X_i = Xr_corr(:,:,pose_i);
            X_j = Xr_corr(:,:,pose_j);
            [err, Ji, Jj] = poseErrorAndJacobian(X_i, X_j, z);
            chi_poses=err'*err;

			chi_stats_poses(i)+=chi_poses;

            pose_i_H=poseMatrixIndex(pose_i, Nr, Nl);
            pose_j_H=poseMatrixIndex(pose_j, Nr, Nl);

            Hr(pose_i_H:pose_i_H+2, pose_i_H:pose_i_H+2)+=Ji'*Ji;
            Hr(pose_i_H:pose_i_H+2, pose_j_H:pose_j_H+2)+=Ji'*Jj;
            Hr(pose_j_H:pose_j_H+2, pose_i_H:pose_i_H+2)+=Jj'*Ji;
            Hr(pose_j_H:pose_j_H+2, pose_j_H:pose_j_H+2)+=Jj'*Jj;

            br(pose_i_H:pose_i_H+2)+=Ji'*err;
            br(pose_j_H:pose_j_H+2)+=Jj'*err;
        end

        Hl = zeros(state_dim, state_dim);
        bl = zeros(state_dim,1);
        #handle landmarks
        for m=1:Ml
            pose_i = associations_Zl_ig(1,m);
            land_j = associations_Zl_ig(2,m);
            z = Zl_ig(1,m);
            Xr_i = Xr_corr(:,:,pose_i);
            Xl_j = Xl_corr(:,land_j);
            [err, Jr_i, Jl_j] = bearingErrorAndJacobian(Xr_i, Xl_j, z);
            chi_bearings=err'*err;

			chi_stats_bearings(i)+=chi_bearings;

            pose_i_H=poseMatrixIndex(pose_i, Nr, Nl);
            land_j_H=landmarkMatrixIndex(land_j, Nr, Nl);

            Hl(pose_i_H:pose_i_H+2, pose_i_H:pose_i_H+2)+=Jr_i'*Jr_i;
            Hl(pose_i_H:pose_i_H+2, land_j_H:land_j_H+1)+=Jr_i'*Jl_j;
            Hl(land_j_H:land_j_H+1, land_j_H:land_j_H+1)+=Jl_j'*Jl_j;
            Hl(land_j_H:land_j_H+1, pose_i_H:pose_i_H+2)+=Jl_j'*Jr_i;

            bl(pose_i_H:pose_i_H+2)+=Jr_i'*err;
            bl(land_j_H:land_j_H+1)+=Jl_j'*err;
        end

        H=Hr+Hl+eye(state_dim,state_dim)*1e-5;
        b=br+bl;

        H((Nr-1)*3+1:Nr*3,:)=[];
        H(:,(Nr-1)*3+1:Nr*3)=[];
        b((Nr-1)*3+1:Nr*3)=[];

        dx = zeros(state_dim,1);
        dx=H\(-b);
        dx = [dx(1:(Nr-1)*3); zeros(3,1); dx((Nr-1)*3+1:end)];
        [Xr_corr, Xl_corr] = box_plus(Xr_corr, Xl_corr, dx, Nr, Nl);
    endfor   
    disp("Done! Optimization completed.");
endfunction