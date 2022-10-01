import math as m
import numpy as nu
import matplotlib.pyplot as plt

pi = m.pi

#### INITIAL SYSTEM PARAMETERS ##################

# geometric parameters of nitrous oxide tank
tank_r = 0.08  # 0.08 m
tank_H = 1.5  # 1.5 m
dip_Tube_H = 0.05  # m
tank_Vol = pi * pow(tank_r, 2) * tank_H

C_d = 0.314
N = 33
A_inj = 0.00000177 * N

A_throat = 0.001  # m
A_exh = 0.00578

# geometric parameters of chamber
init_diameter_Port = 0.040  # m port diameter fuel grain
outer_Diameter_FG = 0.142875  # m outer diameter fuel grain
length_FG = 0.55  # m
a = 0.155
n = 0.555
rho_fuel = 930  # kg/m^3
gas_constant = 385
k_constant = 1.27
combustion_efficiency = 0.95


# Initial conditions
init_Temp = 25  # degC
outlet_Pres = 101300  # Pa

##### SATURATION DATA ###########################
# nitrous saturation pressures(deg C->kPa)
temperatureBreakpoints = [-90.82, -90, -88.46, -85, -80, -75, -70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20,
                          -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 36.42]
pressure = [87.73, 92.29, 101.325, 124.2, 164.2, 213.6, 273.6, 345.7, 431.5, 532.3, 649.9, 785.8, 941.7, 1119, 1321,
            1547, 1801, 2083, 2397, 2744, 3127, 3547, 4007, 4510, 5060, 5660, 6315, 7033, 7251]

# nitrous density table(deg C->Kg / m ^ 3)
liqDensity = [1222.8, 1220.6, 1216.3, 1206.7, 1192.7, 1178.3, 1163.7, 1148.8, 1133.6, 1118.0, 1112.0, 1085.6, 1068.8,
              1051.4, 1033.4, 1014.8, 995.4, 975.2, 953.9, 931.4, 907.4, 881.6, 853.5, 822.2, 786.6, 743.9, 688.0,
              589.4, 452]
vapDensity = [2.613, 2.738, 2.987, 3.609, 4.680, 5.982, 7.546, 9.406, 11.60, 14.16, 17.14, 20.58, 24.53, 29.05, 34.22,
              40.11, 46.82, 54.47, 63.21, 73.26, 84.86, 98.41, 114.5, 133.9, 158.1, 190.0, 236.7, 330.4, 452]

# nitrous latent heat of vapourization(deg C->KJ / Kg)
delHVap = [377, 377, 375, 371, 365, 359, 353, 346, 340, 333, 326, 318, 310, 302, 294, 285, 276, 266, 255, 244, 232, 219,
           204, 188, 169, 147, 117, 64.9, 0]

# specific heat capacity of liquid nitrous(deg C->KJ / Kg * K)
# last two values are sus
specHeatCap = [1.747, 1.750, 1.756, 1.768, 1.781, 1.791, 1.798, 1.803, 1.807, 1.812, 1.818, 1.827, 1.840, 1.858, 1.833,
               1.915, 1.957, 2.011, 2.079, 2.166, 2.274, 2.412, 2.592, 2.834, 3.188, 3.781, 5.143, 100, 1000]

# critical parameters of nox
n2o_pCrit = 7251000;  # critical pressure of n2o[Pa]
n2o_rhoCrit = 452.0;  # critical density of n2o[kg / m ^ 3]
n2o_tCrit = 36.42;  # critical temperature of n2o[deg C]
n2o_ZCrit = 0.28;  # critical compress.factor of n2o[dimless]
n2o_gamma = 1.3;  # gamma of n2o[dimless]

N = 5000
m_v = nu.zeros(N + 1)
m_vc = nu.zeros(N + 1)
del_T = nu.zeros(N + 1)
T_nox = nu.zeros(N + 1)
P_nox = nu.zeros(N + 1)
rho_liq = nu.zeros(N + 1)
rho_vap = nu.zeros(N + 1)
m_liq = nu.zeros(N + 1)
m_vap = nu.zeros(N + 1)
H_v = nu.zeros(N + 1)
c_p_liq = nu.zeros(N + 1)
del_Q = nu.zeros(N + 1)
m_total = nu.zeros(N + 1)
m_dot_oxi = nu.zeros(N + 1)
m_liq_old = nu.zeros(N + 1)
m_liq_new = nu.zeros(N + 1)
m_dot_fuel = nu.zeros(N + 1)
delta_outflow_mass = nu.zeros(N + 1)
m_vap_new = nu.zeros(N + 1)
lagged_m_v = nu.zeros(N + 1)
P_chamber = nu.zeros(N + 1)
OF_ratio = nu.zeros(N + 1)
m_dot_total = nu.zeros(N + 1)
T_chamber = nu.zeros(N + 1)
c_star = nu.zeros(N + 1)
v_exh = nu.zeros(N + 1)
Thrust = nu.zeros(N + 1)
spec_vol_diff = nu.zeros(N + 1)

r_dot = nu.zeros(N + 1)
d_port = nu.zeros(N + 1)
m_fuel_total = nu.zeros(N + 1)

####### MAIN BLOWDOWN CALC #############################

m_v[0] = 0.001  # initial guess
T_nox[0] = init_Temp
P_nox[0] = 1000 * nu.interp(T_nox[0], temperatureBreakpoints, pressure, left=None, right=None, period=None)
rho_liq[0] = nu.interp(T_nox[0], temperatureBreakpoints, liqDensity, left=None, right=None, period=None)
rho_vap[0] = nu.interp(T_nox[0], temperatureBreakpoints, vapDensity, left=None, right=None, period=None)
H_v[0] = nu.interp(T_nox[0], temperatureBreakpoints, delHVap, left=None, right=None, period=None)
c_p_liq[0] = nu.interp(T_nox[0], temperatureBreakpoints, specHeatCap, left=None, right=None, period=None)

init_Liq_Vol = tank_Vol - (pi * pow(tank_r, 2) * dip_Tube_H)
init_Vap_Vol = tank_Vol - init_Liq_Vol
m_liq_init = init_Liq_Vol * rho_liq[0]  # test values
m_total[0] = init_Liq_Vol * rho_liq[0] + init_Vap_Vol * rho_vap[0]
m_liq_new[0] = m_liq_init
m_vap[0] = init_Vap_Vol * rho_vap[0]
P_chamber[0] = 2500000
m_dot_oxi[0] = (C_d * A_inj * pow((2 * rho_liq[0] * (P_nox[0] - P_chamber[0])), 0.5))

init_fuel_mass = pi * length_FG * (pow(outer_Diameter_FG/2, 2) - pow(init_diameter_Port/2, 2)) * rho_fuel
d_port[0] = init_diameter_Port
m_fuel_total[0] = init_fuel_mass
c_star[0] = 1600

lagged_m_v[0] = 0
m_liq_old[0] = m_liq_init

time = 0
del_t = 0.01
i = 0

while m_total[i] > 0.01:
    if m_liq_new[i] > 0.01:
        c_p_liq[i] = nu.interp(T_nox[i], temperatureBreakpoints, specHeatCap, left=None, right=None, period=None)
        H_v[i] = nu.interp(T_nox[i], temperatureBreakpoints, delHVap, left=None, right=None, period=None)

        del_Q[i] = m_v[i] * H_v[i]
        del_T[i + 1] = -(del_Q[i] / (m_liq_old[i] * c_p_liq[i]))
        T_nox[i + 1] = (T_nox[i] + del_T[i])

        P_nox[i + 1] = 1000 * nu.interp(T_nox[i + 1], temperatureBreakpoints, pressure, left=None, right=None, period=None)
        rho_liq[i + 1] = nu.interp(T_nox[i + 1], temperatureBreakpoints, liqDensity, left=None, right=None, period=None)
        rho_vap[i + 1] = nu.interp(T_nox[i + 1], temperatureBreakpoints, vapDensity, left=None, right=None, period=None)

        P_chamber[i + 1] = 2500000
        m_dot_oxi[i + 1] = (C_d * A_inj * pow((2 * rho_liq[i + 1] * (P_nox[i + 1] - P_chamber[i + 1])), 0.5))

        delta_outflow_mass[i + 1] = 0.5 * del_t * (3.0 * m_dot_oxi[i + 1] - m_dot_oxi[i])

        m_total[i + 1] = m_total[i] - delta_outflow_mass[i + 1]
        m_liq_old[i + 1] = m_liq_old[i] - delta_outflow_mass[i + 1]
        spec_vol_diff = ((1 / rho_liq[i + 1]) - (1 / rho_vap[i + 1]))
        m_liq_new[i + 1] = (tank_Vol - (m_total[i + 1] / rho_vap[i + 1])) / spec_vol_diff
        m_vap_new[i + 1] = m_total[i + 1] - m_liq_new[i + 1]
        #####   LAG PROCESS FOR VAPOURIZED MASS    #######
        unlagged_m_v = m_liq_old[i + 1] - m_liq_new[i + 1]
        m_v[i + 1] = m_liq_old[i + 1] - m_liq_new[i + 1]
        tc = del_t / 0.15
        m_v[i + 1] = tc * (m_v[i + 1] - m_liq_old[i + 1]) + m_liq_new[i + 1]
        lagged_m_v[i + 1] = tc * (unlagged_m_v - lagged_m_v[i]) + lagged_m_v[i]
        ##################################################
        m_v[i + 1] = lagged_m_v[i + 1]

        if m_liq_new[i + 1] > m_liq_old[i + 1]:
            m_liq_new[i + 1] = 0
            #break  #ONLY NEEDED IF DOING LIQUID ONLY

        m_liq_old[i + 1] = m_liq_new[i + 1]
    else:
        m_dot_oxi[i] = 1 ##TO BE CONTINUED WITH REST OF VAPOUR BLOWDOWN PORTION (JUST A PLACEHOLDER ATM)

    G_ox = m_dot_oxi[i] / (pi * pow(d_port[i] / 2, 2))
    r_dot[i] = 0.001 * (a * pow(G_ox, n))
    d_port[i + 1] = d_port[i] + 2 * (r_dot[i] * del_t)
    m_dot_fuel[i] = rho_fuel * pi * 2 * length_FG * d_port[i] / 2 * r_dot[i]
    m_fuel_total[i + 1] = m_fuel_total[i] - (del_t * m_dot_fuel[i])
    if d_port[i + 1] > outer_Diameter_FG:
        d_port[i + 1] = outer_Diameter_FG

    OF_ratio[i] = m_dot_oxi[i] / m_dot_fuel[i]
    m_dot_total[i] = m_dot_oxi[i] + m_dot_fuel[i]
    T_chamber[i] = -5742.5 + (7623.5 * OF_ratio[i]) - (2334.4 * pow(OF_ratio[i], 2)) + (314.21 * pow(OF_ratio[i], 3)) \
                   - (15.854 * pow(OF_ratio[i], 4))
    P_chamber[i + 1] = (c_star[i] * m_dot_total[i]) / A_throat
    Mach = pow((T_chamber[i] * gas_constant * k_constant), 0.5)
    c_star[i + 1] = (Mach * combustion_efficiency) / (k_constant * pow((2 / (k_constant + 1)), (k_constant + 1) /
                                                                       (2 * (k_constant - 1))))
    v_exh[i] = pow((((2 * T_chamber[i] * gas_constant * k_constant) / (k_constant - 1)) * (1 - pow((outlet_Pres /
                                                                                                    P_chamber[i]),
                                                                                                   (k_constant - 1)
                                                                                                   / k_constant))), 0.5)
    Mach_exit = v_exh[i] / Mach
    Thrust[i] = (m_dot_total[i] * v_exh[i]) + (outlet_Pres - outlet_Pres) * A_exh

    print('time', time)
    print('T_nox', i, T_nox[i])
    print('P_nox', i, P_nox[i])
    #print('rho_liq', i, rho_liq[i])
    #print('rho_vap', i, rho_vap[i])
    print('m_dot_oxi', i, m_dot_oxi[i])
    print('m_total', i, m_total[i])
    print('m_liq', i, m_liq[i])
    #print('r_dot', i, r_dot[i])
    print('m_dot_fuel', i, m_dot_fuel[i])
    print('d_port', i, d_port[i])
    print('m_fuel_total', i, m_fuel_total[i])
    print('P_chamber', i, P_chamber[i])
    print('c_star', i, c_star[i])
    print('v_exh', i, v_exh[i])
    print('Thrust', i, Thrust[i])
    #print('m_liq_new', i, m_liq_new[i])

    time = time + del_t
    i = i + 1

#while m_vap_new[i] > 0.01:
#    m_vap_new[i + 1] = m_vap_new[i] - 0.01
#    time = time + del_t
#    i = i + 1

plt.plot(Thrust)
plt.ylabel('Thrust')
plt.show()
plt.plot(c_p_liq)
plt.ylabel('c_p_liq')
plt.show()
plt.plot(m_liq_new)
plt.ylabel('m_liq_new')
plt.show()
plt.plot(m_fuel_total)
plt.ylabel('m_fuel_total')
plt.show()
plt.plot(m_total)
plt.ylabel('m_total')
plt.show()
plt.plot(P_nox)
plt.ylabel('P_nox')
plt.show()
plt.plot(T_nox)
plt.ylabel('T_nox')
plt.show()
plt.plot(m_dot_oxi)
plt.ylabel('m_dot_oxi')
plt.show()
plt.plot(m_dot_fuel)
plt.ylabel('m_dot_fuel')
plt.show()
plt.plot(P_chamber)
plt.ylabel('P_chamber')
plt.show()
plt.plot(T_chamber)
plt.ylabel('T_chamber')
plt.show()