# Rate constants
k0 = 0.0039  # min⁻¹
k1 = 0.03    # min⁻¹ a.u.⁻¹
k2 = 0.05    # min⁻¹
k3 = 2.32    # min⁻¹
k4 = 31.14   # min⁻¹ a.u.⁻¹
k5 = 0.32    # min⁻¹
k6 = 2.82    # min⁻¹
k7 = 5.49    # min⁻¹
k8 = 4.37    # min⁻¹
k9 = 18.75   # min⁻¹
k10 = 0.0083 # min⁻¹
k11 = 0.24   # min⁻¹
k12 = 0.30   # min⁻¹
k13 = 1.33   # a.u. min⁻¹
k14 = 0.0048 # a.u.⁻¹
k15 = 0.019  # a.u. min⁻¹
k16 = 0.025  # min⁻¹
k17 = 0.022  # min⁻¹
ke1 = 10     # min⁻¹
ke2 = 0.5    # min⁻¹
kp = ke2*k14/(ke1-k14)
ki = 1/k14
ki1 = 1
ki2 = ki1*ki
k18 = 0.0001

# Scaling factors
WB_scale1 = 3.4
WB_scale2 = 6.6
PCR_scale1 = 0.0026

# Initial conditions
I_tot = 0.0038  # a.u.
x1_0 = 100      # ng/ml
x2_0 = 0        # ng/ml
x3_0 = 0.005    # a.u.
x4_0 = 0        # ng/ml
x5_0 = 0        # ng/ml
x6_0 = 0.5      # ng/ml
x7_0 = 73.6     # a.u.
x8_0 = 0        # a.u.
x9_0 = 1        # ng/ml
x10_0 = I_tot-x2_0  # ng/ml
x11_0 = (x9_0-(x6_0+3*x3_0+x4_0+x5_0))/3  # ng/ml

# Time delays
τ1 = 39   # min
τ2 = 22   # min
τ3 = 164  # min
τ4 = 185  # min

# Distribution parameters
p = 4
q1 = p/τ1
q2 = p/τ2
q3 = p/τ3
q4 = p/τ4 