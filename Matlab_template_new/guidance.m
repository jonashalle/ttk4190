function [chi_d, stop] = guidance(x, y,wp)
persistent k;      % active waypoint index (initialized by: clear LOSchi)
persistent xk yk;  % active waypoint (xk, yk) corresponding to integer k
persistent yp_int;
kappa = 1;
R_switch = 320;
if isempty(k)
    k=1;
    xk = wp(1,k);
    yk = wp(2,k);
    yp_int=0;
end

delta = 1600;
Kp = 1/delta;
Ki = kappa*Kp;
n = length(wp(1,:));

if k<n
    xk_next = wp(1,1+k);
    yk_next = wp(2,1+k);
else
    xk_next = wp(1,n);
    yk_next = wp(2,n);
end

pi_p = atan2((yk_next-yk), (xk_next-xk));
along_track =  (x-xk) * cos(pi_p) + (y-yk) * sin(pi_p);
cross_track =  -(x-xk) * sin(pi_p) + (y-yk) * cos(pi_p);
yp_int_dot = delta*cross_track/(delta^2+(cross_track+kappa*yp_int)^2);

yp_int=yp_int+yp_int_dot;    
distance_wp = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );
stop=false;
if((distance_wp-along_track<R_switch)&&(k<n))
    k = k + 1;
    xk = xk_next;      
    yk = yk_next; 
    
elseif (distance_wp-along_track<100)&&(k==n)
    stop = true;
end

chi_d = ssa(pi_p-atan(Kp*cross_track+Ki*yp_int));