% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
h  = 0.1;    % sampling time [s]
Ns = 10000;    % no. of samples

psi_ref = 10 * pi/180;  
U_ref   = 10;           
psi_ref2 = -20 *pi/180; 
% initial states
eta = [0 0 0]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
x = [nu' eta' delta n]';
L_oa = 161;
A_Lw = 10*L_oa;
V_w = 10;
b_Vw = 135;

a_max = 0.5*pi/180;
r_max = pi/180;
%% PID
w_b = 0.06;
damping = 1;
w_n = 1/(sqrt(1-2*damping^2+sqrt(4*damping^4-4*damping^2+2)))*w_b;
K = 0.0075;
T = 169.5493;
m1 = T/K;
kp = m1*w_n^2;
kd = 2*damping*w_n*m1-1/K;
ki = w_n/10*kp;

x_d= [0 0 0]';
x_ddot = [0 0 0]';
wref = 0.03;
Ad = [0 1 0; 0 0 1; -wref^3 -(2*damping+1)*wref^2 -(2*damping+1)*wref]
Bd = [0 0 wref^3]';
Cd = [1 1 0];
psi_int = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,16);       % table of simulation data
tau_wind =[0 0 0]';
yaw_last = 0;
yaw_int = 0;
for i=1:Ns+1
    
    t = (i-1) * h;              % time (s
    
    
    % current disturbance
    eta_n = x(4:6);
    uc = cosd(45);
    vc = sind(45);
    nu_c = [ uc vc 0 ]';
    %nu_c = [0 0 0]';
%     if t > 400 
%         psi_ref = psi_ref2;
%     end


     
    
    nu_b = Rzyx(0, 0, x(6))'*nu_c;  
    if t > 200 
        % wind disturbance
        u_w = V_w *cosd(b_Vw-x(6));
        v_w = V_w *sind(b_Vw-x(6));
        u_rw = x(1) - u_w;
        v_rw = x(2) - v_w;
        V_rw1 = sqrt(u_rw^2 + v_rw^2); 
        y_rw = -atan2((v_rw),(u_rw));
        C_Y = 0.95*sin(y_rw);
        C_N = 0.15*sin(2*y_rw);
        q = 1/2*1.247*V_rw1^2;
        Ywind = C_Y*A_Lw;
        Nwind = C_N*A_Lw*L_oa;
        tau_wind = q*[0 Ywind Nwind]';
    end
     nu_n = Rzyx(0, 0, x(6))*x(1:3); 

     chi = atan2(nu_n(2), nu_n(1));
     crab_angle = chi - eta_n(3);
     nu_r = x(1:3) - nu_b;
     beta = atan2(nu_r(2), nu_r(1));
   
    % reference models

    x_ddot = Ad*x_d+ Bd*psi_ref;

    psi_d = x_d(1);
    r_d = x_d(2);
    x_d= euler2(x_ddot,x_d,h);
    u_d = U_ref;
    psi_c = psi_d - x(6);
    r_c = r_d -x(3);

    % control law
    psi_int = psi_int + h*psi_c;
    delta_c = (kp*psi_c + ki*psi_int + kd*r_c);
    n_c = 10;             
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c,tau_wind,u_d);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d crab_angle beta];     
 
    % Euler integration
 
    x = euler2(xdot,x,h);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
u       = simdata(:,2);                 % m/s
v       = simdata(:,3);                 % m/s
r       = (180/pi) * simdata(:,4);      % deg/s
x       = simdata(:,5);                 % m
y       = simdata(:,6);                 % m
psi     = (180/pi) * simdata(:,7);      % deg
delta   = (180/pi) * simdata(:,8);      % deg
n       = 60 * simdata(:,9);            % rpm
delta_c = (180/pi) * simdata(:,10);     % deg
n_c     = 60 * simdata(:,11);           % rpm
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
crab_angle   = (180/pi) * simdata(:,15);     % deg
sideslip = (180/pi) * simdata(:,16);     % deg

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');


figure(3)
hold on
plot(t,crab_angle,t,sideslip,'linewidth',2);
hold off
legend("Crab angle","sideslip")
%    title('Crab angle vs. sideslip (deg) (V_c = 0)'); xlabel('time (s)');
title('Crab angle vs sideslip with velocity'); xlabel('time (s)');
grid on