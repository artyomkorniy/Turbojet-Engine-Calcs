clear; clc; close all;

%% Inputs
T_ambient_K        = 288.15;
P_ambient_Pa       = 101325;
compressor_PR      = 1.70;       % estimated, no compressor map
eta_compressor     = 0.72;
combustor_dP_frac  = 0.05;
T_turbine_in_K     = 1150 + 273.15;
eta_turbine        = 0.72;
eta_shaft          = 0.95;
eta_combustion     = 0.98;
fuel_LHV_J_per_kg  = 46.4e6;
mdot_air_kg_s      = 0.15;       % assumed airflow
cp_air_J_kgK       = 1005;
cp_gas_J_kgK       = 1150;
gamma_air          = 1.40;
gamma_gas          = 1.33;
R_air_J_kgK        = 287.05;

assert(T_ambient_K > 0 && P_ambient_Pa > 0 && T_turbine_in_K > 0);
assert(compressor_PR >= 1 && mdot_air_kg_s > 0);
assert(all([eta_compressor,eta_turbine,eta_shaft,eta_combustion] > 0) && ...
       all([eta_compressor,eta_turbine,eta_shaft,eta_combustion] <= 1));
assert(combustor_dP_frac >= 0 && combustor_dP_frac < 1);
assert(cp_air_J_kgK > 0 && cp_gas_J_kgK > 0 && fuel_LHV_J_per_kg > 0);
assert(gamma_air > 1 && gamma_gas > 1 && R_air_J_kgK > 0);

%% Cycle
Tt2 = T_ambient_K;
Pt2 = P_ambient_Pa;
Pt3 = Pt2 * compressor_PR;
Tt3_ideal = Tt2 * compressor_PR^((gamma_air-1)/gamma_air);
Tt3 = Tt2 + (Tt3_ideal-Tt2)/eta_compressor;

Pt4 = Pt3 * (1-combustor_dP_frac);
Tt4 = T_turbine_in_K;

denominator = eta_combustion*fuel_LHV_J_per_kg - cp_gas_J_kgK*Tt4;
assert(denominator > 0, 'Combustor energy balance has no valid positive denominator.');
f = (cp_gas_J_kgK*Tt4 - cp_air_J_kgK*Tt3)/denominator;
assert(f > 0, 'Nonpositive fuel/air ratio: check turbine inlet temperature and parameters.');

specific_compressor_work_J_kg_air = cp_air_J_kgK*(Tt3-Tt2);
Tt5 = Tt4 - specific_compressor_work_J_kg_air/(eta_shaft*(1+f)*cp_gas_J_kgK);
assert(Tt5 > 0 && Tt5 < Tt4, 'Invalid turbine exit temperature.');

Tt5s = Tt4 - (Tt4-Tt5)/eta_turbine;
assert(Tt5s > 0, 'Invalid turbine isentropic exit temperature.');
Pt5 = Pt4*(Tt5s/Tt4)^(gamma_gas/(gamma_gas-1));
assert(Pt5 > 0, 'Invalid turbine exit pressure.');

mdot_fuel_kg_s = f*mdot_air_kg_s;
mdot_exhaust_kg_s = (1+f)*mdot_air_kg_s;
compressor_power_W = mdot_air_kg_s*specific_compressor_work_J_kg_air;
turbine_gas_power_W = mdot_exhaust_kg_s*cp_gas_J_kgK*(Tt4-Tt5);
shaft_power_W = eta_shaft*turbine_gas_power_W;
fuel_chemical_power_W = mdot_fuel_kg_s*fuel_LHV_J_per_kg;

%% Results
fprintf('\nOPEN BRAYTON GAS-PATH MODEL (illustrative inputs)\n');
fprintf('Station       Tt [K]       Pt [kPa]\n');
fprintf('2 (inlet)    %8.2f     %9.2f\n', Tt2,Pt2/1000);
fprintf('3 (comp out) %8.2f     %9.2f\n', Tt3,Pt3/1000);
fprintf('4 (turb in)  %8.2f     %9.2f\n', Tt4,Pt4/1000);
fprintf('5 (turb out) %8.2f     %9.2f\n', Tt5,Pt5/1000);
fprintf('\nFuel/air mass ratio f          = %.5f\n',f);
fprintf('Fuel flow                     = %.5f kg/s\n',mdot_fuel_kg_s);
fprintf('Exhaust flow                  = %.5f kg/s\n',mdot_exhaust_kg_s);
fprintf('Compressor shaft demand       = %.2f kW\n',compressor_power_W/1000);
fprintf('Turbine gas work rate         = %.2f kW\n',turbine_gas_power_W/1000);
fprintf('Useful turbine shaft power    = %.2f kW\n',shaft_power_W/1000);
fprintf('Fuel chemical energy rate     = %.2f kW\n',fuel_chemical_power_W/1000);
if Pt5 <= P_ambient_Pa
    warning('Turbine exit total pressure is <= ambient: a simple expanding nozzle cannot produce the assumed jet.');
end
fprintf('\nNOTE: Compressor PR and air flow are assumptions, not compressor-map results.\n');
fprintf('Turbine inlet temperature is user-supplied; verify whether measured or design target.\n');

%% Plots
stations = {'2 Inlet','3 Compressor','4 Combustor','5 Turbine'};
figure('Name','Turbojet Brayton Gas-Path States','Color','w');
subplot(2,1,1);
plot(1:4,[Tt2,Tt3,Tt4,Tt5],'-o','LineWidth',1.5);
set(gca,'XTick',1:4,'XTickLabel',stations); grid on;
ylabel('Stagnation temperature [K]'); title('Open Brayton Gas-Path States');
subplot(2,1,2);
plot(1:4,[Pt2,Pt3,Pt4,Pt5]/1000,'-o','LineWidth',1.5);
set(gca,'XTick',1:4,'XTickLabel',stations); grid on;
ylabel('Stagnation pressure [kPa]');
