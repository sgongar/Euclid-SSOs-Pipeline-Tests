from numpy import poly1d
from pandas import read_csv

# Reads SSOs catalogue
ssos_df = read_csv('ssos_2.csv', index_col=0)

print(ssos_df.columns)

# # Creates polynomical functions for limits
# upr_coefs_bright = [5.730178e-03, -6.528091e-01, 2.970536e+01,
#                     -6.748832e+02,  7.655361e+03, -3.468232e+04]
# upr_limit_bright = poly1d(upr_coefs_bright)
# lwr_coefs_bright = [5.735456e-03, -6.536302e-01, 2.975250e+01,
#                     -6.761721e+02, 7.672408e+03, -3.477074e+04]
# lwr_limit_bright = poly1d(lwr_coefs_bright)
#
# upr_coefs_faint = [-1.769922e-02, 1.871337e+00, -7.400582e+01,
#                    1.297378e+03, -8.505660e+03]
# upr_limit_faint = poly1d(upr_coefs_faint)
# lwr_coefs_faint = [-3.273145e-02, 3.398693e+00, -1.321923e+02,
#                    2.282417e+03, -1.475840e+04 ]
# lwr_limit_faint = poly1d(lwr_coefs_faint)

# Creates polynomical functions for limits
upr_coefs_bright = [5.730702e-03, -6.528906e-01, 2.971004e+01,
                    -6.750112e+02,  7.657054e+03, -3.469110e+04]
upr_limit_bright = poly1d(upr_coefs_bright)
lwr_coefs_bright = [5.734932e-03, -6.535486e-01, 2.974782e+01,
                    -6.760441e+02, 7.670715e+03, -3.476196e+04]
lwr_limit_bright = poly1d(lwr_coefs_bright)

upr_coefs_faint = [-1.769922e-02, 1.871337e+00, -7.400582e+01,
                   1.297378e+03, -8.505660e+03]
upr_limit_faint = poly1d(upr_coefs_faint)
lwr_coefs_faint = [-3.273145e-02, 3.398693e+00, -1.321923e+02,
                   2.282417e+03, -1.475840e+04 ]
lwr_limit_faint = poly1d(lwr_coefs_faint)



in_bright = 0
out_bright = 0
in_faint = 0
out_faint = 0
# Iterate over sources
for i, row in enumerate(ssos_df.itertuples(), 1):

    if row.MAG_AUTO < 24.5:
        b_low = lwr_limit_bright(float(row.MAG_AUTO))
        b_upr = upr_limit_bright(float(row.MAG_AUTO))

        if b_low < float(row.B_IMAGE) < b_upr:
            in_bright += 1
        else:
            out_bright += 1
    elif row.MAG_AUTO > 24.5:
        b_low = lwr_limit_faint(float(row.MAG_AUTO))
        b_upr = upr_limit_faint(float(row.MAG_AUTO))

        if b_low < float(row.B_IMAGE) < b_upr:
            in_faint += 1
        else:
	    out_faint += 1

print('bright -- in {} - out {}'.format(in_bright, out_bright))
print('faint  -- in {} - out {}'.format(in_faint, out_faint))

