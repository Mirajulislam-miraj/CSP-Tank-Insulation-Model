clc; clear; close all;
 
%% Geometry (one tank)
D = 32;
r_inner = D/2;
H = 15;
 
A_side   = 2*pi*r_inner*H;
A_top    = pi*r_inner^2;
A_bottom = A_top;
 
%% Temperatures (K)
T_hot  = 565 + 273.15;
T_cold = 290 + 273.15;
T_amb  = 283;
 
%% Molten salt NANO3 60% & KNO3 40%
m_salt = 2.4106e7;
Cp     = 1500;
 
E_tank = m_salt * Cp * (T_hot - T_cold);
 
fprintf('Total Thermal capacity = %.0f MWh\n', E_tank/3.6e9);
 
%% Insulation stack [thickness (m), k (W/m-K)]
layers = [
     0.035  45;     % Carbon steel shell (35 mm)
    0.40   0.35;   % IFB
    0.15   0.04    % Mineral wool
];
 
%% Radial conduction resistance (SIDE)
r_outer = r_inner;
R_side = 0;
 
for j = 1:size(layers,1)
    t = layers(j,1);
    k = layers(j,2);
    r_in = r_outer;
    r_outer = r_outer + t;
    R_side = R_side + log(r_outer/r_in)/(2*pi*k*H);
end
 
%% Planar conduction resistance (TOP & BOTTOM)
R_top = 0;
for j = 1:size(layers,1)
    t = layers(j,1);
    k = layers(j,2);
    R_top = R_top + t/(k*A_top);
end
R_bottom = R_top;
 
%% External heat transfer
h_conv  = 8;
epsilon = 0.8;
sigma   = 5.67e-8;
 
%% Time
dt = 3600;
Nt = 24;
time_h = (0:Nt-1)';
 
%% Initial state (FIXED)
T = T_hot;       % Start at hot temperature
Q = E_tank;      % Stored energy at hot condition
 
Temp = zeros(Nt,1);
Energy = zeros(Nt,1);
Eff = zeros(Nt,1);
Loss_conv = zeros(Nt,1);
Loss_rad  = zeros(Nt,1);
 
%% ===================== TIME LOOP =====================
for i = 1:Nt
 
    %% ---- SIDE WALL (solve Ts via Newton–Raphson)
    Ts = T - 20;
 
    for it = 1:30
        Q_cond_side = (T - Ts)/R_side;
        Q_conv_side = h_conv*A_side*(Ts - T_amb);
        Q_rad_side  = epsilon*sigma*A_side*(Ts^4 - T_amb^4);
 
        f  = Q_cond_side - (Q_conv_side + Q_rad_side);
        df = -1/R_side - h_conv*A_side ...
             - 4*epsilon*sigma*A_side*Ts^3;
 
        Ts_new = Ts - f/df;
        if abs(Ts_new - Ts) < 1e-6
            Ts = Ts_new;
            break;
        end
        Ts = Ts_new;
    end
 
    Q_side = Q_conv_side + Q_rad_side;
 
    %% ---- TOP SURFACE (solve Ts via Newton–Raphson)
    Ts_top = T - 20;
 
    for it = 1:30
        Q_cond_top = (T - Ts_top)/R_top;
        Q_conv_top = h_conv*A_top*(Ts_top - T_amb);
        Q_rad_top  = epsilon*sigma*A_top*(Ts_top^4 - T_amb^4);
 
        f  = Q_cond_top - (Q_conv_top + Q_rad_top);
        df = -1/R_top - h_conv*A_top ...
             - 4*epsilon*sigma*A_top*Ts_top^3;
 
        Ts_new = Ts_top - f/df;
        if abs(Ts_new - Ts_top) < 1e-6
            Ts_top = Ts_new;
            break;
        end
        Ts_top = Ts_new;
    end
 
    Q_top = Q_conv_top + Q_rad_top;
 
    %% ---- BOTTOM SURFACE (solve Ts via Newton–Raphson)
    Ts_bot = T - 20;
 
    for it = 1:30
        Q_cond_bot = (T - Ts_bot)/R_bottom;
        Q_conv_bot = h_conv*A_bottom*(Ts_bot - T_amb);
        Q_rad_bot  = epsilon*sigma*A_bottom*(Ts_bot^4 - T_amb^4);
 
        f  = Q_cond_bot - (Q_conv_bot + Q_rad_bot);
        df = -1/R_bottom - h_conv*A_bottom ...
             - 4*epsilon*sigma*A_bottom*Ts_bot^3;
 
        Ts_new = Ts_bot - f/df;
        if abs(Ts_new - Ts_bot) < 1e-6
            Ts_bot = Ts_new;
            break;
        end
        Ts_bot = Ts_new;
    end
 
    Q_bot = Q_conv_bot + Q_rad_bot;
 
    %% ---- TOTAL HEAT LOSS
    Q_loss = (Q_side + Q_top + Q_bot) * dt;
 
    %% ---- UPDATE ENERGY
    Q = max(Q - Q_loss, 0);
 
    %% ---- UPDATE TEMPERATURE (FIXED)
    T = T_cold + Q/(m_salt*Cp);
    T = max(T, T_cold);
 
    %% ---- STORE
    Temp(i)   = T - 273.15;
    Energy(i) = Q/3.6e9;
    Eff(i)    = 100 * Q/E_tank;
 
    Loss_conv(i) = ((Q_conv_side + Q_conv_top + Q_conv_bot)*dt)/3.6e6;
    Loss_rad(i)  = ((Q_rad_side  + Q_rad_top  + Q_rad_bot )*dt)/3.6e6;
end
 
%% ===================== SUMMARY =====================
Loss_1tank = sum(Loss_conv + Loss_rad)/1000;
 
fprintf('\n24 h results:\n');
fprintf('Hot tank loss  = %.3f MWh\n', Loss_1tank);
fprintf('Hot Tank efficiency = %.2f %%\n', Eff(end));
 
subplot(2,2,1)
plot(time_h,Temp,'b','LineWidth',2)
grid on; xlabel('Time (h)'); ylabel('Temperature (°C)')
title('Tank Temperature')
 
subplot(2,2,2)
plot(time_h,Energy,'r','LineWidth',2)
grid on; xlabel('Time (h)'); ylabel('Energy (MWh)')
title('Stored Energy (Hot Tank)')
 
subplot(2,2,3)
plot(time_h,Loss_conv,'g','LineWidth',2); hold on
plot(time_h,Loss_rad,'m','LineWidth',2)
plot(time_h,Loss_conv+Loss_rad,'k--','LineWidth',2)
grid on; xlabel('Time (h)'); ylabel('Loss (kWh)')
legend('Convection','Radiation','Total')
title('Heat Loss Breakdown (Hot Tank)')
 
subplot(2,2,4)
plot(time_h,Eff,'c','LineWidth',2)
grid on; xlabel('Time (h)'); ylabel('Efficiency (%)')
title('Efficiency (Hot Tank) ')