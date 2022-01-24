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
