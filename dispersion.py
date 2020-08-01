from cmath import sqrt, nan

# omega=np.zeros([101,6])
i = 0
j = 0
omega = [[0 for i in range(101)] for i in range(6)]

for j in range(0, 6):

    if (j == 0):
        n = 2
        print("Equatorial inertia-gravity wave 2nd merdional mode")
        for i in range(0, 101):
            k = -5 + i * 0.1
            omega[j][i] = sqrt(k * k + 2 * n + 1.0)


    elif (j == 1):
        n = 1
        print("Equatorial inertia-gravity wave 1st merdional mode")
        for i in range(0, 101):
            k = -5 + i * 0.1
            omega[j][i] = sqrt(k * k + 2 * n + 1.0)

    elif (j == 2):
        n = 0
        print("Equatorial mixed Rossby-gravity wave")
        for i in range(0, 101):
            k = -5 + i * 0.1
            omega[j][i] = k / 2.0 + sqrt((k / 2.0) * (k / 2.0) + 1.0)

    elif (j == 3):
        n = 1
        print("Equatorial Rossby wave 1st merdional mode")
        for i in range(0, 101):
            k = -5 + i * 0.1
            omega[j][i] = -k / (k * k + 2 * n + 1.0)

    elif (j == 4):
        n = 2
        print("Equatorial Rossby wave 2st merdional mode")
        for i in range(0, 101):
            k = -5 + i * 0.1
            omega[j][i] = -k / (k * k + 2 * n + 1.0)

        for r in range(0, 30):
            for i in range(0, 101):
                k = -5 + i * 0.1
                om = omega[j][i]
                try:
                    dom = -(om * om - (k * k + 2.0 * n + 1.0) - k / om) / (2.0 * om + k / om / om)
                except:
                    dom = nan

                omega[j][i] = omega[j][i] + dom


    elif (j == 5):
        print("Equatorial Kelvin wave")
        for i in range(0, 101):
            n = -1
            k = -5 + i * 0.1
            omega[j][i] = k

for i in range(0, 101):
    k = -5 + i * 0.1
    print('%.10e' % (k), '%.10e' % (omega[0][i].real), '%.10e' % (omega[1][i].real),
          '%.10e' % (omega[2][i].real), '%.10e' % (omega[3][i].real), '%.10e' % (omega[4][i].real),
          '%.10e' % (omega[5][i].real))
    print('\n')


