1;

function A = v2t(v)
	c = cos(v(3));
	s = sin(v(3));

	A=[c, -s, v(1);
	   s,  c, v(2);
	   0,  0,   1];
endfunction;


function v = t2v(A)
	v(1:2,1) = A(1:2,3);
	v(3,1) = atan2(A(2,1), A(1,1));
endfunction;


function v=flattenIsometryByColumns(T)
        v=zeros(12,1);
        v(1:9)=reshape(T(1:3,1:3),9,1);
        v(10:12)=T(1:3,4);
endfunction


function delta_theta = box_minus(theta1, theta2)
	difference = theta2-theta1;
    if (difference>pi)
        delta_theta = abs(difference-2*pi);
	elseif (difference<-pi)
        delta_theta = abs(difference+2*pi);
    else
        delta_theta = abs(difference);
	endif
endfunction;


function point = intersect_lines(p11,p12,p21,p22)
        k1 = (p12(2)-p11(2))/(p12(1)-p11(1));
        k2 = (p22(2)-p21(2))/(p22(1)-p21(1));
        w21 = p11(2)-k1*p11(1);
        w22 = p21(2)-k2*p21(1);
        x=(w22-w21)/(k1-k2);
        y=w21+k1*x;
        point = [x;y];
endfunction;


function [e, Ji, Jj]=poseErrorAndJacobian(Xi,Xj,Z)
        Ri=Xi(1:2,1:2);
        Rj=Xj(1:2,1:2);

        ti=Xi(1:2,3);
        tj=Xj(1:2,3);

        dR_0=[0 -1;1 0];

        Z_hat=eye(3);
        Z_hat(1:2,1:2)=Ri'*Rj;
        Z_hat(1:2,3)=Ri'*(tj-ti);

        e=flattenIsometryByColumns(Z_hat-Z);

        Ji=zeros(6,3);
        Jj=zeros(6,3);

        Jj(1:4,3)=reshape(Ri'*dR_0*Rj, 4, 1);
        Jj(5:6,1:2)=Ri';

        Jj(5:6,3)=Ri'*dR_0*tj;
        Ji=-Jj;
endfunction







function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)

        pose_dim=6;
        landmark_dim=3;

        if (pose_index>num_poses)
                v_idx=-1;
                return;
        endif;

        v_idx=1+(pose_index-1)*pose_dim;
endfunction;

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
       
        pose_dim=6;
        landmark_dim=3;

        if (landmark_index>num_landmarks)
                v_idx=-1;
                return;
        endif;
        
        v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;