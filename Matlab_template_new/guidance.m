function [chi_d] = guidance(x, y,wp)
persistent k;      % active waypoint index (initialized by: clear LOSchi)
persistent xk yk;  % active waypoint (xk, yk) corresponding to integer k

R_switch = 300;
if isempty(k)
    k=1;
    xk = wp(1,k);
    yk = wp(2,k);
    
end

delta = 1000;
Kp = 1/delta;

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

distance_wp = sqrt( (xk_next-xk)^2 + (yk_next-yk)^2 );

if((distance_wp-along_track<R_switch)&&(k<n))
    k = k + 1;
    xk = xk_next;       % update active waypoint
    yk = yk_next; 
end

chi_d = ssa(pi_p-atan(Kp*cross_track));