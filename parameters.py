import numpy as np

# COMMON PARAMETERS
X_MAX_SELF = 0.70
X_MAX_REC = 0.70
X_MAX_YUMI = 0.90
COS_EPS = 0.05

m_arrow = 0.025
eta_d = 0.80

# SELF-BOW
ell_self = 0.75
y0_self = 0.10
theta0_self = 20 * np.pi/180
k_self = 620.0
L0_self = y0_self + ell_self * np.cos(theta0_self)

# RECURVE
ell_rec = 0.65
y0_rec = 0.10
theta0_rec = 25 * np.pi/180
k_rec = 620.0
ell_r = 0.03
alpha = 2.0
L0_rec = y0_rec + ell_rec * np.cos(theta0_rec)

# YUMI
ell1_yumi = 0.40
ell2_yumi = 0.88
y0_yumi = 0.15
theta0_yumi = 20 * np.pi/180
k1_yumi = 650.0
k2_yumi = 590.0

ell_avg_yumi = (ell1_yumi + ell2_yumi) / 2
k_avg_yumi = (k1_yumi + k2_yumi) / 2
L0_yumi = y0_yumi + ell_avg_yumi * np.cos(theta0_yumi)
