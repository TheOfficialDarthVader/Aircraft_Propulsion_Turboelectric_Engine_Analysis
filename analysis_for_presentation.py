# -*- coding: utf-8 -*-

# Import the required libraries
# --------------------------------------------------------------------------- #
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# --------------------------------------------------------------------------- #

# Properties taken from the NEA specifictation
# --------------------------------------------------------------------------- #
c_p_air = 1005  # J/(K kg)
gamma = 1.4
a = 296.508  # m/s
P_a = 23836  # Pa
V_flight = 231.276  # m/s
M_flight = 0.78

P_02 = 35570  # Pa
T_02 = 245.4  # K
T_023 = 327.0  # K
T_03 = 793.5  # K
T_04 = 1500  # K
T_045 = 1033.5  # K
T_05 = 577.7  # K

thrust_required = 38900  # N
lift_over_drag = 21.6
LCV = 43000000  # J
fan_istrp_eff = 0.9
F_N = 38.9 * 1000   # N
V_jc = V_jbp = 340  # m/s
m_dot_fan_inlet = 392  # kg/s
m_dot_core = 35  # kg/s
m_dot_fuel = 0.441/3600/9.81*thrust_required  # kg/s
non_dim_mass_flow_rate = 1.078  # Assumes M=0.6 fan inlet
# --------------------------------------------------------------------------- #

# Variables used for formatting the plots
# --------------------------------------------------------------------------- #
tick_spacing_fans = 1
tick_spacing_bpr = 1
# --------------------------------------------------------------------------- #

# -------------------------- Bypass Ratio Analysis -------------------------- #
# Initialise arrays for storing values to be plotted
# --------------------------------------------------------------------------- #
fpr_vals = np.arange(1.3, 1.8+0.01, 0.01)
bpr_vals = []
bpr_vals_var = []
sfc_vals = []
# --------------------------------------------------------------------------- #
# Loop through the fpr values to calculate the associated bpr and sfc
# --------------------------------------------------------------------------- #
for fpr in fpr_vals:
    # Calculate T_013 assuming ideal gas an using isentropic efficiency
    T_013_s = T_02 * (fpr)**((gamma-1)/(gamma))
    T_013 = T_02 + ((T_013_s-T_02)/(fan_istrp_eff))
    P_013 = P_02 * fpr

    # Calculate the bypass pressure ratio
    bpr = ((m_dot_core*c_p_air*(T_045-T_05)) -
           (m_dot_core*c_p_air*(T_023-T_02))) / \
          (m_dot_core*c_p_air*(T_013-T_02))
    bpr_vals.append(bpr)

    # Calculate the bpyass jet velocity
    V_jbp = np.sqrt(2*c_p_air*T_013*(1-(P_a/P_013)**((gamma-1)/(gamma))))

    # Calculate the specific net thrust
    F_N_specific = (V_jc-V_flight) + bpr*(V_jbp-V_flight)

    # Calcualte the overall efficiency
    overall_eff = (V_flight*F_N_specific/1000)/(c_p_air*(T_04-T_03))

    # Calculate the specific fuel consumption
    sfc_vals.append((V_flight/(overall_eff*LCV))*1000)
# --------------------------------------------------------------------------- #
# Plot the results of bpr and sfc
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
fig.tight_layout(pad=30.0)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(fpr_vals, bpr_vals, color='black')
ax1.set_xlim(min(fpr_vals), max(fpr_vals))
ax1.set_ylim(round(min(bpr_vals))-2, max(bpr_vals)+1)
ax1.set_xlabel('Fan Pressure Ratio')
ax1.set_ylabel('Bypass Ratio')
ax1.set_title('Bypass ratio as a Function of Fan Pressure Ratio')

ax2.plot(bpr_vals, sfc_vals, color='black')
ax2.set_xlim(round(min(bpr_vals)), round(max(bpr_vals)))
ax2.set_ylim(round(min(sfc_vals))-0.2, max(sfc_vals)+0.2)
ax2.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_bpr))
ax2.set_xlabel('Bypass Pressure Ratio')
ax2.set_ylabel(r'Specific Fuel Consumption, $g/(s KN)$')
ax2.set_title('Specific Fuel Consumption as a Function of Bypass Pressure ' +
              'Ratio')
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# ------------------------- Number of fans Analysis ------------------------- #
# Initialise arrays for storing values to be plotted
# --------------------------------------------------------------------------- #
fan_vals = np.arange(1, 20, 0.11)
num_fans = 1
weight_vals = []
sfc_corrected = []
bpr = 11.2
# --------------------------------------------------------------------------- #
# Loop through the fan values to calculate the associated weight and efficiency
# --------------------------------------------------------------------------- #
for num_fans in fan_vals:
    # Calculate the mass flow through 1 fan
    m_dot_fan = ((m_dot_core*bpr))/num_fans

    # Calculate the diameter of the fan required
    d_fan = 2*np.sqrt((m_dot_fan*np.sqrt(c_p_air*T_02)) /
                      (np.pi*non_dim_mass_flow_rate*P_02))

    # Calculate the weight of a single fan
    k = 12/((3)**(2.4))
    weight_fan = k*d_fan**(2.4)*1000

    # Calculate the weight of the propulsion system
    weight_total_fans = weight_fan * num_fans
    weight_vals.append(weight_total_fans)

    # Calculate the bare engine thrust
    F_N_bare = m_dot_fan_inlet*(V_jc-V_flight)

    # Calculate engine specific thrust
    X = V_jc-V_flight

    # Calculate the effective thrust
    F_N_effective = F_N_bare*(1-(9.25/X))

    # Calcualte the corrected thrust
    F_N_corrected = F_N_effective - weight_total_fans/(lift_over_drag)

    # Calcualte the corrected specific fuel consumption
    sfc_corrected.append((1000*1000*m_dot_fuel)/F_N_corrected)
# --------------------------------------------------------------------------- #
# Plot the results of weight and specific fuel consumption
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
fig.tight_layout(pad=30.0)
ax3 = fig.add_subplot(211)
ax4 = fig.add_subplot(212)

ax3.plot(fan_vals, weight_vals, color='black')
ax3.set_xlim(1, max(fan_vals))
ax3.set_ylim(min(weight_vals), max(weight_vals))
ax3.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_fans))
ax3.set_xlabel('Number of Fans')
ax3.set_ylabel(r'Total Propulsion System Weight (1 engine), $kg$')
ax3.set_title('Propulsion System Weight as a Function of Number of Fans')

ax4.plot(fan_vals, sfc_corrected, color='black')
ax4.set_xlim(1, max(fan_vals))
ax4.set_ylim(12.50, 12.6)
ax4.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_fans))
ax4.set_xlabel('Number of Fans')
ax4.set_ylabel(r'Corrected $sfc$, $g/(s KN)$')
ax4.set_title(r'Corrected $sfc$ as a Function of Number of Fans')
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# -------------------------- Nacelle Drag Analysis -------------------------- #
# Initialise arrays for storing values to be plotted
# --------------------------------------------------------------------------- #
nacelle_drag_vals = []
# --------------------------------------------------------------------------- #
# Loop through the fan values to calculate the associated nacelle drag
# --------------------------------------------------------------------------- #
for num_fans in fan_vals:
    k = 0.04
    nacelle_drag_vals.append(k*V_flight*(F_N/X))
# --------------------------------------------------------------------------- #
# Plot the results of nacelle drag
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
fig.tight_layout(pad=30.0)
ax5 = fig.add_subplot(211)

ax5.plot(fan_vals, nacelle_drag_vals, color='black')
ax5.set_xlim(1, max(fan_vals))
ax5.set_ylim(min(nacelle_drag_vals), max(nacelle_drag_vals))
ax5.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_fans))
ax5.set_xlabel('Number of Fans')
ax5.set_ylabel(r'Nacelle drag, $N$')
ax5.set_title('Nacelle Drag as a Function of Number of Fans')
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------- Viability Analysis ---------------------------- #
# Initialise arrays for storing values to be plotted
# --------------------------------------------------------------------------- #
sfc_corrected_super_vals = []
sfc_corrected_non_super_vals = []
# --------------------------------------------------------------------------- #
# Set the power densities of the super and non-superconducting systems
# --------------------------------------------------------------------------- #
power_density_super = 36000  # From NASA paper
power_density_non_super = 6000  # Assumed from assignment instructions
# --------------------------------------------------------------------------- #
# Loop through the fan values to calculate the associated nacelle drag
# --------------------------------------------------------------------------- #
for fan_weight in weight_vals:
    # Calculat the power required for the aircraft to fly at it's cruise speed
    power_required = V_flight*thrust_required

    # Calculate the weights of the systems
    weight_of_super = fan_weight + power_required/power_density_super
    weight_of_non_super = fan_weight + power_required/power_density_non_super

    # Calculate the net force from the bare engine
    F_N_bare = m_dot_fan_inlet*(V_jc-V_flight)
    X = V_jc-V_flight

    # Calculate the effective force on the engine (takes into account drag)
    F_N_effective = F_N_bare*(1-(9.25/X))

    # Calculate the corrected force (takes into account weight)
    F_N_corrected_super = F_N_effective - (weight_of_super /
                                           lift_over_drag)
    F_N_corrected_non_super = F_N_effective - (weight_of_non_super /
                                               lift_over_drag)

    # Calculate the corrected specific fuel consumption
    sfc_corrected_super_vals.append((1000*1000*m_dot_fuel) /
                                    F_N_corrected_super)
    sfc_corrected_non_super_vals.append((1000*1000*m_dot_fuel) /
                                        F_N_corrected_non_super)
# --------------------------------------------------------------------------- #
# Plot the results of corrected sfc
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
fig.tight_layout(pad=30.0)
ax6 = fig.add_subplot(211)

ax6.plot(fan_vals, sfc_corrected_super_vals, color='blue')
ax6.plot(fan_vals, sfc_corrected_non_super_vals, color='green')
ax6.set_xlim(1, max(fan_vals))
ax6.set_ylim(12.52,
             max(sfc_corrected_non_super_vals)+0.01)
ax6.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_fans))
ax6.set_xlabel('Number of Fans')
ax6.set_ylabel(r'Corrected $sfc$, $g/(s KN)$')
ax6.set_title(r'Corrected $sfc$ as a Function of Number of Fans')
ax6.legend(['Super conducting electronics', 'Non-superconducting electronics'])
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
