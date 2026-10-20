clear; clc; close all;

%% Inputs
T_ambient_K       = 288.15;
P_ambient_Pa      = 101325;
V_flight_m_s      = 0.0;
compressor_PR     = 1.70;    % estimated
eta_compressor    = 0.72;
combustor_dP_frac = 0.05;
T_turbine_in_K    = 1150 + 273.15;
eta_turbine       = 0.72;
eta_shaft         = 0.95;
eta_combustion    = 0.98;
fuel_LHV_J_kg     = 46.4e6;
cp_air_J_kgK      = 1005;
cp_gas_J_kgK      = 1150;
gamma_air         = 1.40;
gamma_gas         = 1.33;
R_gas_J_kgK       = cp_gas_J_kgK*(gamma_gas-1)/gamma_gas;
nozzle_exit_D_m   = 3*0.0254; % assumes exit is the throat
assumed_compressor_airflow_kg_s = 0.15;  % not measured
Cd_nozzle         = 0.98;
Cv_nozzle         = 0.98;

assert(T_ambient_K>0 && P_ambient_Pa>0 && T_turbine_in_K>0);
assert(V_flight_m_s>=0 && compressor_PR>=1 && nozzle_exit_D_m>0);
assert(assumed_compressor_airflow_kg_s>0);
assert(all([eta_compressor,eta_turbine,eta_shaft,eta_combustion,Cd_nozzle,Cv_nozzle]>0) && ...
       all([eta_compressor,eta_turbine,eta_shaft,eta_combustion,Cd_nozzle,Cv_nozzle]<=1));
assert(combustor_dP_frac>=0 && combustor_dP_frac<1);
assert(cp_air_J_kgK>0 && cp_gas_J_kgK>0 && gamma_air>1 && gamma_gas>1 && R_gas_J_kgK>0);
assert(fuel_LHV_J_kg>0);

%% Cycle
assert(V_flight_m_s==0, 'This stationary-test-stand model requires V_flight_m_s = 0.');
Tt2 = T_ambient_K; Pt2 = P_ambient_Pa;
Pt3 = Pt2*compressor_PR;
Tt3s = Tt2*compressor_PR^((gamma_air-1)/gamma_air);
Tt3 = Tt2+(Tt3s-Tt2)/eta_compressor;
Pt4 = Pt3*(1-combustor_dP_frac); Tt4 = T_turbine_in_K;
denom = eta_combustion*fuel_LHV_J_kg-cp_gas_J_kgK*Tt4;
assert(denom>0,'No valid positive denominator in fuel energy balance.');
f = (cp_gas_J_kgK*Tt4-cp_air_J_kgK*Tt3)/denom;
assert(f>0, 'Invalid fuel/air ratio; check input parameters.');
Tt5 = Tt4-cp_air_J_kgK*(Tt3-Tt2)/(eta_shaft*(1+f)*cp_gas_J_kgK);
assert(Tt5>0 && Tt5<Tt4,'Invalid turbine work balance.');
Tt5s = Tt4-(Tt4-Tt5)/eta_turbine;
assert(Tt5s>0,'Invalid turbine pressure expansion.');
Pt5 = Pt4*(Tt5s/Tt4)^(gamma_gas/(gamma_gas-1));

%% Nozzle
assert(Pt5>P_ambient_Pa, 'Nozzle total pressure must exceed ambient to expand.');
Ae = pi*(nozzle_exit_D_m^2)/4;
critical_static_to_total = (2/(gamma_gas+1))^(gamma_gas/(gamma_gas-1));
back_to_total = P_ambient_Pa/Pt5;
is_choked = back_to_total <= critical_static_to_total;
if is_choked
    Pe = Pt5*critical_static_to_total;
    Te = Tt5*2/(gamma_gas+1);
else
    Pe = P_ambient_Pa;
    Te = Tt5*(Pe/Pt5)^((gamma_gas-1)/gamma_gas);
end

Ve_ideal = sqrt(max(0, 2*cp_gas_J_kgK*(Tt5-Te)));
Ve = Cv_nozzle*Ve_ideal;

mdot_exhaust = Cd_nozzle*Ae*(Pe/(R_gas_J_kgK*Te))*Ve_ideal;
assert(mdot_exhaust>0,'No positive exhaust mass flow.');
mdot_air = mdot_exhaust/(1+f);
mdot_fuel = mdot_exhaust-mdot_air;

F_momentum_N = mdot_exhaust*Ve-mdot_air*V_flight_m_s;
F_pressure_N = (Pe-P_ambient_Pa)*Ae;
F_net_N = F_momentum_N+F_pressure_N;

%% Results
fprintf('\nTHEORETICAL STATIC TURBOJET THRUST: ILLUSTRATIVE INPUTS\n');
fprintf('Compressor pressure ratio   = %.3f\n',compressor_PR);
fprintf('Combustor exit Tt4          = %.1f K\n',Tt4);
fprintf('Turbine exit Pt5            = %.2f kPa\n',Pt5/1000);
fprintf('Turbine exit Tt5            = %.2f K\n',Tt5);
fprintf('Fuel/air mass ratio         = %.5f\n',f);
fprintf('Convergent nozzle exit D    = %.2f mm\n',1000*nozzle_exit_D_m);
fprintf('Nozzle choked?              = %d (1=yes, 0=no)\n',is_choked);
fprintf('Exit static pressure        = %.2f kPa\n',Pe/1000);
fprintf('Exit static temperature     = %.2f K\n',Te);
fprintf('Exit velocity               = %.2f m/s\n',Ve);
fprintf('Air mass flow (NOZZLE-implied) = %.5f kg/s\n',mdot_air);
fprintf('Illustrative compressor flow = %.5f kg/s\n',assumed_compressor_airflow_kg_s);
fprintf('Nozzle/compressor flow ratio = %.2f\n',mdot_air/assumed_compressor_airflow_kg_s);
fprintf('Fuel mass flow (INFERRED)   = %.5f kg/s\n',mdot_fuel);
fprintf('Momentum contribution       = %.2f N\n',F_momentum_N);
fprintf('Pressure contribution       = %.2f N\n',F_pressure_N);
fprintf('NOZZLE-BASED CONDITIONAL THRUST = %.2f N\n',F_net_N);
if abs(mdot_air/assumed_compressor_airflow_kg_s - 1) > 0.2
    warning(['Nozzle-based flow differs by >20% from the illustrative compressor flow. ', ...
             'The assumed pressure ratio, flow and 3-inch exit are NOT matched; ', ...
             'do NOT treat the output as a physically consistent turbojet thrust estimate.']);
end
fprintf('\nWARNING: This is NOT measured thrust; inputs are assumed.\n');
fprintf('No compressor/turbine maps were available; actual operating point is UNKNOWN.\n');
fprintf('The nozzle-implied flow is NOT a measurement or a validated turbocharger prediction.\n');
fprintf('The 3-inch diameter must be the nozzle MINIMUM area for this nozzle model.\n');
