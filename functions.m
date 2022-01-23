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
endfunction

function delta_theta = box_minus(theta1, theta2)
	difference = theta2-theta1;
    if (difference>pi)
        delta_theta = abs(difference-2*pi);
	elseif (difference<-pi)
        delta_theta = abs(difference+2*pi);
    else
        delta_theta = abs(difference);
endfunction