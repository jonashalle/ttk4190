% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;    % sampling time [s]
Ns = 60000;    % no. of samples
load WP.mat
psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_ref   = 9;            % desired surge speed (m/s)

%Constants
rho_a = 1.247;
Vc = 1;
Bv_c = 45;
Bv_w = 135;
Vw = 10;
Loa = 161;
ALw = 10*Loa;

% initial states
eta = [0 0 -110*pi/180]';  %Posisjon
nu  = [0.1 0 0]'; %Hastighet
delta = 0;
n = 0;
Qm = 0;
x = [nu' eta' delta n Qm]';





% PID controller
omega_b = 0.06;
zeta = 1;
T = 169.546;
K = 0.0075;
m = T/K;
d = 1/K;
omega_n = 1/(sqrt(1-2*zeta^2+sqrt(4*zeta^4-4*zeta^2+2)))*omega_b;
kp = m*omega_n^2;
ki = omega_n/10*kp;
kd = 2*zeta*omega_n*m-d;
pid_params = [kp ki kd];

%reference model
xd = [0 0 0]';
xd_dot = [0 0 0]';
omega_ref = 0.03;
Ad = [0 1 0;
      0 0 1;
      -omega_ref^3 -(2*zeta+1)*omega_ref^2 -(2*zeta+1)*omega_ref];
Bd = [0 0 omega_ref^3]';
Cd = [1 1 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,18);       % table of simulation data
tau_wind = [0 0 0]';
psi_integral = 0;
for i=1:Ns+1
    
    t = (i-1) * h;              % time (s)
    
    % current disturbance
    eta_n = x(4:6);
    uc = cosd(Bv_c)*Vc;
    vc = sind(Bv_c)*Vc;
    nu_c = [ uc vc 0 ]';
    
    % wind disturbance
    if t > 200   

        uw = Vw*cosd(Bv_w-x(6));
        vw = Vw*sind(Bv_w-x(6));
        u_rw = x(1)-uw;
        v_rw = x(2)-vw;
        Vrw = sqrt(u_rw^2+v_rw^2);
    
        gamma_rw = -atan2(v_rw,u_rw);
        C_Y = 0.95*sin(gamma_rw);
        C_N = 0.15*sin(2*gamma_rw);
        
        Ywind = C_Y*ALw;
        Nwind = C_N*ALw*Loa;
        
        tau_wind = 0.5*rho_a*Vrw^2*[0 Ywind Nwind]';
    end
    %Heading and crab
     nu_b = Rzyx(0, 0, x(6))'*nu_c;  
     nu_n = Rzyx(0, 0, x(6))*x(1:3); 
     chi = atan2(nu_n(2), nu_n(1));
     crab_angle = chi - eta_n(3);
     nu_r = x(1:3) - nu_b;
     beta = atan2(nu_r(2), nu_r(1));
    
%     if t > 400
%         psi_ref = -20 * pi/180;  % desired yaw angle (rad)
%     end
    [chi_d,stop] = guidance(x(4),x(5),WP);
    if stop == true
        U_ref = 0;
    end
    psi_ref = chi_d;
    %psi_ref = chi_d-crab_angle;
    
    % reference models
    xd_dot = Ad*xd+Bd*psi_ref;
    xd = euler2(xd_dot,xd,h);
    psi_d = xd(1);
    r_d = xd(2);
    u_d = U_ref;
    
    psi_c = psi_d - x(6);
    r_c = r_d-x(3);
    % control law
    psi_integral = psi_integral + psi_c*h;
    wind_up_lim = 5;
    if psi_integral<-wind_up_lim
        psi_integral=-wind_up_lim;

    elseif psi_integral>wind_up_lim
        psi_integral=wind_up_lim;
    end 
    n_c = 10;                % propeller speed (rps)
     
    delta_c = (kp*psi_c+kd*r_c+ki*psi_integral);              % rudder angle command (rad)

    % ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c,tau_wind,U_ref);
    
    % store simulation data in a table (for testing)    
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d crab_angle beta chi_d chi];   
    
    % Euler integration
    x = euler2(xdot,x,h);    
    %add reference
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
chi_d = (180/pi) * simdata(:,17);     % deg
chi = (180/pi) * simdata(:,18);     % deg

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

figure(4)
hold on
plot(t,chi,t,chi_d,'linewidth',2)
hold off;
legend('chi','chi desired')
title('Actual and commanded course'); xlabel('time (s)');
pathplotter(x,y)
hold off