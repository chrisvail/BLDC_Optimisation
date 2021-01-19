clc
clear all

%% Defining Constants
global RHO_LI;
global RHO_C;
global RHO_M;
global RHO_MS;
global RHO_A;
RHO_LI = 8050;
RHO_C = 8940;
RHO_M = 7007;
RHO_MS = 7850;
RHO_A = 2710;

global I_max;
global V;
global cruise_torque;
global cruise_omega;
global torque_multiple;
global magnet_strength;
I_max = 5;
V = 12;
cruise_torque = 0.04;
cruise_omega = 11000*pi/30;
torque_multiple = 10;
magnet_strength = 1;

%% Optimisation Set Up
% Getting objective function

% Input variables
% 1   n_windings      integer            
% 2   n_teeth         integer           
% 3   tooth_width     continuous               
% 4   height          continuous               
% 5   r_stator        continuous          
% 6   r_wire          continuous          
% 7   t_air_gap       continuous        
% 8   t_magnets       continuous  

% Bounds
lower_bounds = [5, 3, 2E-03, 10E-03, 10E-03, 0.08E-03, 0.4E-03, 2E-03];
upper_bounds = [100, 15, 10E-03, 50E-03, 50E-03, 5E-03, 0.4E-03, 2E-03];

% Linear constraints
A=[];
Aeq = [];
beq = [];
b = [];

objective_function = @get_motor_mass;

%% Optimisation
% Genetic Algorithm
disp("GA:");
options = optimoptions(@ga, 'Display', 'final', 'ConstraintTolerance', 1E-07);%, 'PlotFcn',{@gaplotbestf,@gaplotstopping}); % 
[x_ga, fval_ga] = ga(objective_function, 8, A, b, Aeq, beq, lower_bounds, upper_bounds, @nonlinear_constraints, [1, 2], options);



x0 = x_ga;
% SQP
disp("SQP:");
options = optimoptions(@fmincon, 'Algorithm','sqp','MaxFunctionEvaluations', 1000);
[x_sqp, fval_sqp] = fmincon(objective_function,x0,A,b,Aeq,beq,lower_bounds,upper_bounds, @nonlinear_constraints, options);


% Global Search
disp("Global Search");
gs = GlobalSearch;
options = optimoptions(@fmincon, 'Algorithm','sqp');
problem = createOptimProblem("fmincon", "objective", objective_function, "x0", x0, "lb", lower_bounds, "ub", upper_bounds, "nonlcon", @nonlinear_constraints,"options", options);
[x_gs, f_gs] = run(gs, problem);


disp("x* - GA");
disp([1 1 1000 1000 1000 1000 1000 1000].*x_ga);
disp(fval_ga*1000 + "g");
disp("x* - SQP");
disp([1 1 1000 1000 1000 1000 1000 1000].*x_sqp);
disp(fval_sqp*1000 + "g");
disp("x* - GS");
disp([1 1 1000 1000 1000 1000 1000 1000].*x_gs);
disp(f_gs*1000 + "g");


%% Objective Function
function mass = get_motor_mass(x)
    global RHO_LI;
    global RHO_C;
    global RHO_M;
    global RHO_MS;
    global RHO_A;
    
    tooth_weight = (x(3)*x(4)*(x(5)*0.8 - 2.5E-3) + 2.5E-3*x(4)*((2*pi*x(5))/(3*x(2)) - (2*pi*x(5)/200)))*RHO_LI + hollow_cylinder_weight(0.1*x(5), 0.2*x(5), x(4), RHO_LI);
    winding_weight = cylinder_weight(x(6), x(1)*(2*x(3) + 2*x(4)), RHO_C);
    magnet_weight = hollow_cylinder_weight(x(5)+x(7), x(5)+x(7) + x(8), x(4), RHO_M);
    rotor_weight = hollow_cylinder_weight(x(5)+x(7) + x(8), x(5)+x(7) + x(8)*2, x(4), RHO_MS);
    can_weight = 2*cylinder_weight(x(5) + x(7) + (x(8)*2), 1E-03, RHO_A);
    
    mass = tooth_weight + winding_weight + magnet_weight + rotor_weight + can_weight;
end

function mass = get_motor_mass_vector(x)
    global RHO_LI;
    global RHO_C;
    global RHO_M;
    global RHO_MS;
    global RHO_A;
    
    tooth_weight = (x(3)*x(4)*(x(5)*0.8 - 2.5E-3) + 2.5E-3*x(4)*((2*pi*x(5))/(3*x(2)) - (2*pi*x(5)/200)))*RHO_LI  + hollow_cylinder_weight(0.1*x(5), 0.2*x(5), x(4), RHO_LI);
    winding_weight = cylinder_weight(x(6), x(1)*(2*x(3) + 2*x(4)), RHO_C);
    magnet_weight = hollow_cylinder_weight(x(5)+x(7), x(5)+x(7) + x(8), x(4), RHO_M);
    rotor_weight = hollow_cylinder_weight(x(5)+x(7) + x(8), x(5)+x(7) + x(8)*2, x(4), RHO_MS);
    can_weight = 2*cylinder_weight(x(5)+x(7) + x(8)*2, 4E-03, RHO_A);
    
    mass = [tooth_weight winding_weight magnet_weight rotor_weight can_weight];
end

function w = cylinder_weight(r, h, rho)
    w = pi*r*r*h*rho;
end

function w = hollow_cylinder_weight(r_i, r_o, h, rho)
    w = pi*h*rho*(r_o^2 - r_i^2);
end

%% Non-linear Constraints

function [c, ceq] = nonlinear_constraints(x)
    c(1) = tooth_width_constraint(x);
    c(2) = winding_constraint(x);
    c(3) = power_constraint(x);
    c(4) = torque_constraint(x);
    ceq = [];
end

% Tooth width constraint
function constraint = tooth_width_constraint(x)
    total_tooth_width = 3*x(2)*x(3);
    r_0 = x(5)*0.4;
    circumference = 2*pi*r_0;
    constraint = circumference - total_tooth_width;
end

% Winding Constraint
function constraint = winding_constraint(x)
    r_stator = x(5) - 2.5E-3;
    r_wire = x(6);
    r_0 = r_stator*0.4;
    
    theta = pi/(3*x(2));
    
    total_windings = 0;
    for i=0:2
        winding_layer = (r_stator - max((x(3)/2 + r_wire + i*sqrt(3)*r_wire)/tan(theta), r_0))/(2*r_wire);
        if winding_layer > 0
            total_windings = total_windings + winding_layer;
        end
    end
    constraint = x(1) - total_windings;
end

% Power Constraint
function constraint = power_constraint(x)
    constraint = get_min_operating_current(x) - get_max_operating_current(x);
end

function current = get_min_operating_current(x)
    global V;
    global cruise_torque;
    global cruise_omega;
    R = get_winding_resistance(x);
    discriminant = V^2 - 4*R*cruise_torque*cruise_omega;
    if discriminant < 0
        % set infeasible current
        current = 10000;
    else
        current = (V - sqrt(discriminant))/(2*R);
    end
end

function R = get_winding_resistance(x)
    winding_length = get_winding_length(x);
    winding_area = pi*x(6)^2;
    R = winding_length*1.68E-6/winding_area;
end

function winding_length = get_winding_length(x)
    winding_length = 2*x(1)*x(2)*(2*x(3) + 2*x(4));
end

function current = get_max_operating_current(x)
    global I_max;
    % Empirically generated model
    max_wire_current = 8.73E06*x(6)^2 + 281*x(6) - 4.28E-02;
    current = min(max_wire_current, I_max);
end

% Torque Constraint
function constraint = torque_constraint(x)
    constraint = get_min_required_torque() - get_max_produced_torque(x);
end

function torque = get_min_required_torque()
    global cruise_torque;
    global torque_multiple;
    torque = cruise_torque*torque_multiple;
end

function torque = get_max_produced_torque(x)
    global magnet_strength;
    B = magnet_strength*x(8)/(x(7) + x(8));
    max_current = get_max_operating_current(x);
    winding_length = get_winding_length(x);
    
    torque = 4*B*max_current*winding_length*x(5);
end