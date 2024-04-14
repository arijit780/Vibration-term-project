# Vibration-term-project

TUNED MASS DAMPER MATLAB CODE :

%% Plots

plot_mode_shapes = true;

%% Building parameters

mass_floor = 2e6;
EI_column = 41325000;
Height_column = 3.3;
num_floors = 20;
    
SA_all = true; % Add or remove shock absorbers at every floor
SA = false; % Add or remove shock absorber at the bottom. Note that SA_all adds one at the bottom. Use this if only one at the bottom is needed
c_sa = 4e5*4;

tmd = true; % Turn tmd on or off
tmd_floor = 20;
m_tmd = 2e6*5;
c_tmd = 4e6;
k_tmd = 1e6;

earthquake_amp = 0.5; %m
earthquake_freq = 2; %Hz


%% Load definition

time_interval = 0.01;
time_vector = 0:time_interval:1000;
num_steps = 1000/time_interval;
num_steps_wind_ramp_up = 30/time_interval;
num_steps_wind_ramp_down = num_steps_wind_ramp_up + 30/time_interval;
num_steps_earthquake = 10/time_interval;

wind_loading = zeros(1,num_steps+1);
for i = 1:num_steps+1
    if i < num_steps_wind_ramp_up+1
        wind_loading(i) = time_vector(i)/30;
    elseif i < num_steps_wind_ramp_down+1
        wind_loading(i) = -time_vector(i)/30 + 2;
    else
        wind_loading(i) = 0;
    end
end
clearvars i
wind_floor_multiplier = 375.*(1:num_floors)';
wind_loading = wind_loading.*wind_floor_multiplier;
clearvars wind_floor_multiplier

earthquake_loading = zeros(1,num_steps+1);
for i = 1:num_steps+1
    if i < num_steps_earthquake+1
        earthquake_loading(i) = (24*EI_column/Height_column^3)*earthquake_amp*cos(earthquake_freq*2*pi*time_vector(i)) + SA*c_sa*(0.5)*(-earthquake_amp*earthquake_freq*2*pi*sin(earthquake_freq*2*pi*time_vector(i)));
    end
end
clearvars i
earthquake_loading = [earthquake_loading;zeros(num_floors-1,num_steps+1)];

total_force = wind_loading + earthquake_loading;

if tmd
    total_force = [total_force; zeros(1,size(time_vector,2))];
end

%% Matrices

Mass_matrix = mass_floor.*eye(num_floors);
if tmd
    Mass_matrix(end+1,end+1) = m_tmd;
end

K_matrix = 2.*eye(num_floors);
for i = 1:num_floors-1
    K_matrix(i,i+1) = -1;
    K_matrix(i+1,i) = -1;
end
clearvars i;
K_matrix(end,end) = 1;
K_matrix = (24*EI_column/Height_column^3).*K_matrix;
if tmd
    K_matrix(tmd_floor, tmd_floor) = K_matrix(tmd_floor, tmd_floor) + k_tmd;
    K_matrix(end+1, tmd_floor) = -k_tmd;
    K_matrix(tmd_floor, end + 1) = -k_tmd;
    K_matrix(end,end) = k_tmd;
end

C_matrix = zeros(num_floors);
if SA_all
    C_matrix = eye(num_floors);
    for i = 1:num_floors-1
        C_matrix(i,i+1) = -0.5;
        C_matrix(i+1,i) = -0.5;
    end
    clearvars i;
    C_matrix(end,end) = 0.5;
    C_matrix = c_sa.*C_matrix;
end
if tmd
    C_matrix(tmd_floor, tmd_floor) = C_matrix(tmd_floor, tmd_floor) + c_tmd;
    C_matrix(end+1, tmd_floor) = -c_tmd;
    C_matrix(tmd_floor, end + 1) = -c_tmd;
    C_matrix(end,end) = c_tmd;
end
if SA
    C_matrix(1,1) = c_sa/2;
end


%% Natural modes and frequencies

[V, Lambda] = eig(K_matrix, Mass_matrix);
U = V.*sqrt(mass_floor);
Omega = sqrt(Lambda);
nat_freq = diag(Omega);
nat_freq_Hz = nat_freq./(2*pi);

if(plot_mode_shapes)
    figure(1);
    for i = 1:20
        subplot(4,5,i);
        plot([0;U(1:num_floors,i)], 0:num_floors, 'b-', 'Linewidth', 1.5);
        title("Mode " + num2str(i) + ", fn = " + num2str(nat_freq_Hz(i)) + "Hz");
        grid on;
        xlim([-0.8;0.8]);
        hold on;
    end
    clearvars i;
end


%% Response
A = [zeros(num_floors + tmd), eye(num_floors+tmd); -Mass_matrix\K_matrix, -Mass_matrix\C_matrix];
B = [zeros(num_floors+tmd); inv(Mass_matrix)];

Phi = exp(1)^(A*time_interval);
Gamma = A\(Phi - eye(2*(num_floors+tmd)))*B;

resp = zeros(2*(num_floors+tmd), size(time_vector,2));

for i = 1:num_steps
    resp(:,i+1) = Phi*resp(:,i) + Gamma*total_force(:,i); 
end

plot(time_vector, resp(20,:), 'LineWidth', 2)
hold on;
plot(time_vector, resp(20 + num_floors,:), 'LineWidth', 2) % Assuming the response for the TMD is stored right after the main structure response
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Response with and without TMD')
legend('Without TMD', 'With TMD')
grid on;
xlim([0 100]);




IC Engine matlab code:


%% Engine parameters
bore = 0.08; % Bore (m)
stroke = 0.09; % Stroke (m)
conn_rod_len = 0.14; % Connecting rod length (m)
rpm = 2000; % Engine speed (rpm)
mass_rec = 0.5; % Reciprocating mass (kg)
mass_rot = 2.5; % Rotating mass (kg)
radius = stroke / 2; % Crank radius (m)
n_cycles = 10; % Number of cycles to simulate

%% Flywheel parameters
flywheel_inertia = 0.05; % Flywheel moment of inertia (kg-m^2)
eccentric_mass = 0.02; % Eccentric mass (kg)
eccentric_distance = 0.01; % Distance of eccentric mass from center (m)

%% Simulation parameters
t_step = 0.0001; % Time step (s)
theta_step = 2 * pi / 1000; % Crank angle step (rad)
time = 0:t_step:(n_cycles * 2 * pi / (rpm / 60)); % Time array

%% Compute engine constants
V_disp = pi * bore^2 * stroke / 4; % Displaced volume (m^3)
r_cr = conn_rod_len / radius; % Connecting rod ratio
omega = rpm * 2 * pi / 60; % Angular velocity (rad/s)

%% Initialize arrays
theta = 0:theta_step:(2 * pi * n_cycles); % Crank angle array
speed_no_flywheel = zeros(size(time)); % Speed without flywheel
speed_with_flywheel = zeros(size(time)); % Speed with flywheel
unbalanced_force = zeros(size(theta)); % Unbalanced force

%% Compute unbalanced force and speed fluctuation
for i = 1:length(time)
    theta_idx = mod(floor(time(i) * omega / theta_step), length(theta)) + 1;
    theta_curr = theta(theta_idx);
    sin_theta = sin(theta_curr);
    cos_theta = cos(theta_curr);

% Compute unbalanced force due to reciprocating mass and eccentric mass
    unbalanced_force_rec = mass_rec * radius * omega^2 * (cos_theta + cos_theta / r_cr);
    unbalanced_force_ecc = eccentric_mass * eccentric_distance * omega^2 * cos(theta_curr);
    unbalanced_force(theta_idx) = unbalanced_force_rec + unbalanced_force_ecc;
   
  
%% Engine parameters
bore = 0.08; % Bore (m)
stroke = 0.09; % Stroke (m)
conn_rod_len = 0.14; % Connecting rod length (m)
rpm = 2000; % Engine speed (rpm)
mass_rec = 0.5; % Reciprocating mass (kg)
mass_rot = 2.5; % Rotating mass (kg)
radius = stroke / 2; % Crank radius (m)
n_cycles = 10; % Number of cycles to simulate

%% Flywheel parameters
flywheel_inertia = 0.05; % Flywheel moment of inertia (kg-m^2)
eccentric_mass = 0.02; % Eccentric mass (kg)
eccentric_distance = 0.01; % Distance of eccentric mass from center (m)

%% Simulation parameters
t_step = 0.0001; % Time step (s)
theta_step = 2 * pi / 1000; % Crank angle step (rad)
time = 0:t_step:(n_cycles * 2 * pi / (rpm / 60)); % Time array

%% Compute engine constants
V_disp = pi * bore^2 * stroke / 4; % Displaced volume (m^3)
r_cr = conn_rod_len / radius; % Connecting rod ratio
omega = rpm * 2 * pi / 60; % Angular velocity (rad/s)

%% Initialize arrays
theta = 0:theta_step:(2 * pi * n_cycles); % Crank angle array
speed_no_flywheel = zeros(size(time)); % Speed without flywheel
speed_with_flywheel = zeros(size(time)); % Speed with flywheel
unbalanced_force = zeros(size(theta)); % Unbalanced force

%% Compute unbalanced force and speed fluctuation
for i = 1:length(time)
    theta_idx = mod(floor(time(i) * omega / theta_step), length(theta)) + 1;
    theta_curr = theta(theta_idx);
    sin_theta = sin(theta_curr);
    cos_theta = cos(theta_curr);

% Compute unbalanced force due to reciprocating mass and eccentric mass
    unbalanced_force_rec = mass_rec * radius * omega^2 * (cos_theta + cos_theta / r_cr);
    unbalanced_force_ecc = eccentric_mass * eccentric_distance * omega^2 * cos(theta_curr);
    unbalanced_force(theta_idx) = unbalanced_force_rec + unbalanced_force_ecc;
   
 
       % Compute speed without flywheel
   speed_no_flywheel(i) = omega * (1 + cos_theta / r_cr);
   
    % Compute speed with flywheel
   angular_accel = unbalanced_force(theta_idx) * radius / (mass_rot * radius^2 + flywheel_inertia);
    if i > 1
        speed_with_flywheel(i) = speed_with_flywheel(i-1) + angular_accel * t_step;
    else
        speed_with_flywheel(i) = omega;
    end
end

%% Plot results
figure;
subplot(2, 1, 1);
plot(time, speed_no_flywheel, time, speed_with_flywheel);
title('Engine Speed');
xlabel('Time (s)');
ylabel('Speed (rad/s)');
legend('Without Flywheel', 'With Flywheel');
subplot(2, 1, 2);
plot(theta, unbalanced_force);
title('Unbalanced Force');
xlabel('Crank Angle (rad)');
ylabel('Force (N)');
end

%% Plot results
figure;
subplot(2, 1, 1);
plot(time, speed_no_flywheel, time, speed_with_flywheel);
title('Engine Speed');
xlabel('Time (s)');
ylabel('Speed (rad/s)');
legend('Without Flywheel', 'With Flywheel');
subplot(2, 1, 2);
plot(theta, unbalanced_force);
title('Unbalanced Force');
xlabel('Crank Angle (rad)');
ylabel('Force (N)');
