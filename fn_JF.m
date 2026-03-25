%% partial_h/partial_X : Jacobian of measurement eq
function F = fn_JF(xhat, rates, dt)
A = zeros(3,3);

phi = xhat(1);
theta = xhat(2);

p = rates(1);
q = rates(2);
r = rates(3);

A(1,1) = q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta);
A(1,2) = q*sin(phi)*sec(theta)^2 + r*cos(phi)*sec(theta)^2;
A(1,3) = 0;

A(2,1) = -q*sin(phi) - r*cos(phi);
A(2,2) = 0; 
A(2,3) = 0;

A(3,1) = q*cos(phi)*sec(theta) - r*sin(phi)*sec(theta);
A(3,2) = q*sin(phi)*sec(theta)*tan(theta) + r*cos(phi)*sec(theta)*tan(theta);
A(3,3) = 0;

F = eye(3,3) + A*dt;