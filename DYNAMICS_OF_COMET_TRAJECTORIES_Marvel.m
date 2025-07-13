% Simulation of Binary Star System with Planet Orbiting Star 2 and Comet

% Constants
G = 6.674e-11; % Gravitational constant (m^3 kg^-1 s^-2)

% Star Parameters (Binary System)
M1 = 1.5e30; % Mass of Star 1 (kg)
M2 = 1.0e30; % Mass of Star 2 (kg)
total_mass = M1 + M2;

% Binary Star System - Center of Mass Adjustment
r1 = [1.0e11, 0];  % Initial position of Star 1 
r2 = [-2.0e11, 0]; % Adjusted position of Star 2 
com = (M1 * r1 + M2 * r2) / total_mass; % Center of mass/ Barycenter
r1 = r1 - com; % Adjust position of Star 1 respect to CoM
r2 = r2 - com; % Adjust position of Star 2 respect to CoM

% Compute stable orbital velocities for both stars
r12 = norm(r1 - r2);
v_orb = sqrt(G * total_mass / r12); % Orbital velocity for binary stars
v1 = [0, (M2 / total_mass) * v_orb]; % Orbital velocity of Star 1
v2 = [0, (-M1 / total_mass) * v_orb]; % Orbital velocity of Star 2

% Planet Parameters (Orbiting Star 2)
Mp = 6.0e24; % Mass of the planet (kg)
r_p = [4.0e11, 0]; % Adjusted distance from Star 2 for a stable orbit
rp = r2 + r_p; % Absolute position of the planet around Star 2

% **Ensure Circular Orbit for Planet**
% Circular velocity formula
r_p_distance = norm(rp - r2); % Orbit radius around Star 2
vp_mag = sqrt(G * M2 / r_p_distance); % Correct orbital velocity for a circular orbit
vp = [0, vp_mag] + v2; % Adjusted velocity for the planet orbiting Star 2

% Simulation Parameters
time_step = 1e4; % Smaller time step (s) for more precision
total_time = 1.0e8; % Total simulation time (s)
n_steps = total_time / time_step; % Number of steps

% Initialize position and velocity arrays
pos_star1 = zeros(n_steps, 2); % 2D positions for Star 1
pos_star2 = zeros(n_steps, 2); % 2D positions for Star 2
pos_planet = zeros(n_steps, 2); % 2D positions for the planet
pos_comet = zeros(n_steps, 2); % 2D positions for the comet

pos_star1(1, :) = r1;
pos_star2(1, :) = r2;
pos_planet(1, :) = rp;

v_star1 = v1;
v_star2 = v2;
v_planet = vp;

% Comet Parameters (Entering the System from a Different Direction)
mc = 3.0e14; % Mass of the comet (kg)

% Starting position of the comet farther and entering from a different direction
rc = [5.0e11, 3.5e11]; % Comet starts at a position farther to the right
% Set comet’s initial velocity vector towards the center of mass (com) with high speed
direction_to_comet = com - rc; % Direction vector pointing to the center of mass
direction_to_comet = direction_to_comet / norm(direction_to_comet); % Normalize
vc = 4.0e4 * direction_to_comet; % Increased speed towards the center of mass (entering fast)

% Assign comet initial position and velocity
pos_comet(1, :) = rc;
v_comet = vc;

% Initialize velocity magnitude arrays
v_star1_mag = zeros(n_steps, 1);
v_star2_mag = zeros(n_steps, 1);
v_planet_mag = zeros(n_steps, 1);
v_comet_mag = zeros(n_steps, 1);

% --- Define collision detection radii ---
collision_radius_star = 1e10;    % Approximate radius for stars
collision_radius_planet = 6.4e6; % Approximate Earth radius

% Initialize collision flag
collision_happened = false;
collision_time = 0;
collision_position = [0, 0];
collision_target = '';

% Simulation Loop with Collision Detection
for i = 2:n_steps
    % Calculate distances between bodies
    d12 = norm(pos_star1(i - 1, :) - pos_star2(i - 1, :));
    d1p = norm(pos_star1(i - 1, :) - pos_planet(i - 1, :));
    d2p = norm(pos_star2(i - 1, :) - pos_planet(i - 1, :));
    d1c = norm(pos_star1(i - 1, :) - pos_comet(i - 1, :));
    d2c = norm(pos_star2(i - 1, :) - pos_comet(i - 1, :));
    dpc = norm(pos_planet(i - 1, :) - pos_comet(i - 1, :));

    % Gravitational forces
    F12 = G * M1 * M2 / d12^3 * (pos_star2(i - 1, :) - pos_star1(i - 1, :));
    F1p = G * M1 * Mp / d1p^3 * (pos_planet(i - 1, :) - pos_star1(i - 1, :));
    F2p = G * M2 * Mp / d2p^3 * (pos_planet(i - 1, :) - pos_star2(i - 1, :));
    F1c = G * M1 * mc / d1c^3 * (pos_comet(i - 1, :) - pos_star1(i - 1, :));
    F2c = G * M2 * mc / d2c^3 * (pos_comet(i - 1, :) - pos_star2(i - 1, :));
    Fpc = G * Mp * mc / dpc^3 * (pos_comet(i - 1, :) - pos_planet(i - 1, :));

    % Update velocities using total forces
    a_star1 = (F12 + F1p + F1c) / M1;
    a_star2 = (-F12 + F2p + F2c) / M2;
    a_planet = (-F1p - F2p + Fpc) / Mp;
    a_comet = (-F1c - F2c - Fpc) / mc;

    v_star1 = v_star1 + a_star1 * time_step;
    v_star2 = v_star2 + a_star2 * time_step;
    v_planet = v_planet + a_planet * time_step;
    v_comet = v_comet + a_comet * time_step;

    % Update positions
    pos_star1(i, :) = pos_star1(i - 1, :) + v_star1 * time_step;
    pos_star2(i, :) = pos_star2(i - 1, :) + v_star2 * time_step;
    pos_planet(i, :) = pos_planet(i - 1, :) + v_planet * time_step;
    pos_comet(i, :) = pos_comet(i - 1, :) + v_comet * time_step;

    % Compute velocity magnitudes
    v_star1_mag(i) = norm(v_star1);
    v_star2_mag(i) = norm(v_star2);
    v_planet_mag(i) = norm(v_planet);
    v_comet_mag(i) = norm(v_comet);

    % --- Collision Detection ---
    if ~collision_happened
        if d1c <= collision_radius_star
            fprintf('Collision: Comet hit Star 1 at time %.2f seconds.\n', (i-1)*time_step);
            collision_happened = true;
            collision_time = (i-1) * time_step;
            collision_position = pos_comet(i, :);
            collision_target = 'Star 1';
        elseif d2c <= collision_radius_star
            fprintf('Collision: Comet hit Star 2 at time %.2f seconds.\n', (i-1)*time_step);
            collision_happened = true;
            collision_time = (i-1) * time_step;
            collision_position = pos_comet(i, :);
            collision_target = 'Star 2';
        elseif dpc <= collision_radius_planet
            fprintf('Collision: Comet hit the Planet at time %.2f seconds.\n', (i-1)*time_step);
            collision_happened = true;
            collision_time = (i-1) * time_step;
            collision_position = pos_comet(i, :);
            collision_target = 'Planet';
        end
    end

    % Stop motion of comet on collision
    if collision_happened
        v_comet = [0, 0];  % Freeze comet
    end
end
% === Velocity Magnitude Subplots ===
figure;

% Plot 1: Velocity graph for Star 1
subplot(2, 2, 1);
plot((1:n_steps)*time_step, v_star1_mag, 'r');
title('Velocity of Star 1');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot 2: Velocity graph for Star 2
subplot(2, 2, 2);
plot((1:n_steps)*time_step, v_star2_mag, 'b');
title('Velocity of Star 2');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot 3: Velocity graph for Planet
subplot(2, 2, 3);
plot((1:n_steps)*time_step, v_planet_mag, 'g');
title('Velocity of Planet');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot 4: Velocity graph for Comet
subplot(2, 2, 4);
plot((1:n_steps)*time_step, v_comet_mag, 'k');
title('Velocity of Comet');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
% === Animation ===
figure;
hold on;
axis equal;
xlim([-4e11, 4e11]);
ylim([-4e11, 4e11]);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Simulation of Binary Star System with Planet and Comet');
grid on;

plot(pos_star1(:, 1), pos_star1(:, 2), 'r:', 'DisplayName', 'Star 1 Path');
plot(pos_star2(:, 1), pos_star2(:, 2), 'b:', 'DisplayName', 'Star 2 Path');
plot(pos_planet(:, 1), pos_planet(:, 2), 'g-', 'DisplayName', 'Planet Path');
plot(pos_comet(:, 1), pos_comet(:, 2), 'k-', 'DisplayName', 'Comet Path');

star1_marker = plot(pos_star1(1, 1), pos_star1(1, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Star 1');
star2_marker = plot(pos_star2(1, 1), pos_star2(1, 2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Star 2');
planet_marker = plot(pos_planet(1, 1), pos_planet(1, 2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Planet');
comet_marker = plot(pos_comet(1, 1), pos_comet(1, 2), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Comet');

legend;

for i = 1:n_steps
    set(star1_marker, 'XData', pos_star1(i, 1), 'YData', pos_star1(i, 2));
    set(star2_marker, 'XData', pos_star2(i, 1), 'YData', pos_star2(i, 2));
    set(planet_marker, 'XData', pos_planet(i, 1), 'YData', pos_planet(i, 2));
    set(comet_marker, 'XData', pos_comet(i, 1), 'YData', pos_comet(i, 2));
    
    % Change comet color after collision
    if collision_happened && i*time_step >= collision_time
        set(comet_marker, 'MarkerFaceColor', 'r');
    end
    
    drawnow;
end

% Display collision marker
if collision_happened
    plot(collision_position(1), collision_position(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Star 1 Path', 'Star 2 Path', 'Planet Path', 'Comet Path', ...
           'Star 1', 'Star 2', 'Planet', 'Comet', 'Collision');
end

hold off;
