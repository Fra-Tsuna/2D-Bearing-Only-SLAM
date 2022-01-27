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

    for i=1:iterations
        H = zeros(state_dim, state_dim);
        b = zeros(state_dim,1);

        %handle poses
        for m=1:Mr
            omega = eye(6);
            omega(1:4,1:4) = (1e-3)*omega(1:4,1:4);
            pose_i = associations_Zr_ig(1,m);
            pose_j = associations_Zr_ig(2,m);
            z = Zr_ig(:,:,m);
            X_i = Xr_corr(:,:,pose_i);
            X_j = Xr_corr(:,:,pose_j);
            [err, Ji, Jj] = poseErrorAndJacobian(X_i, X_j, z);
            chi_poses=err'*omega*err;
			chi_stats_poses(i)+=chi_poses;

            pose_i_H=poseMatrixIndex(pose_i, Nr, Nl);
            pose_j_H=poseMatrixIndex(pose_j, Nr, Nl);

            H(pose_i_H:pose_i_H+3-1,
            pose_i_H:pose_i_H+3-1)+=Ji'*omega*Ji;

            H(pose_i_H:pose_i_H+3-1,
            pose_j_H:pose_j_H+3-1)+=Ji'*omega*Jj;

            H(pose_j_H:pose_j_H+3-1,
            pose_i_H:pose_i_H+3-1)+=Jj'*omega*Ji;

            H(pose_j_H:pose_j_H+3-1,
            pose_j_H:pose_j_H+3-1)+=Jj'*omega*Jj;

            b(pose_i_H:pose_i_H+3-1)+=Ji'*omega*err;
            b(pose_j_H:pose_j_H+3-1)+=Jj'*omega*err;
        end

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

            H(pose_i_H:pose_i_H+3-1,
            pose_i_H:pose_i_H+3-1)+=Jr_i'*Jr_i;

            H(pose_i_H:pose_i_H+3-1,
            land_j_H:land_j_H+2-1)+=Jr_i'*Jl_j;

            H(land_j_H:land_j_H+2-1,
            land_j_H:land_j_H+2-1)+=Jl_j'*Jl_j;

            H(land_j_H:land_j_H+2-1,
            pose_i_H:pose_i_H+3-1)+=Jl_j'*Jr_i;

            b(pose_i_H:pose_i_H+3-1)+=Jr_i'*err;
            b(land_j_H:land_j_H+2-1)+=Jl_j'*err;
        end

        dx = zeros(state_dim,1);
        H((Nr-1)*3+1:Nr*3,:)=[];
        H(:,(Nr-1)*3+1:Nr*3)=[];
        b((Nr-1)*3+1:Nr*3)=[];
        dx=-H\b;
        dx = [dx(1:(Nr-1)*3)' zeros(1,3) dx((Nr-1)*3+1:end)']';
        [Xr_corr, Xl_corr] = box_plus(Xr_corr, Xl_corr, dx, Nr, Nl) 
    endfor   
endfunction