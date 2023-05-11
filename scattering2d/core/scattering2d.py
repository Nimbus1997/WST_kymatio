# modi3
# date: 2023.05.09
# what chagne?
# 3) change to make the S_3 as well
# - changed the 'filter_bank.py' too (only make one psi filter for each 'j' (no 'res' difference) - for simplisty )


def scattering2d(x, pad, unpad, backend, J, L, phi, psi, max_order,
        out_type='array'):
    subsample_fourier = backend.subsample_fourier
    modulus = backend.modulus
    rfft = backend.rfft
    ifft = backend.ifft
    irfft = backend.irfft    
    cdgmm = backend.cdgmm
    concatenate = backend.concatenate

    # Define lists for output.
    # out_S_0, out_S_1, out_S_2 = [], [], []
    out_S_0, out_S_1, out_S_2, out_S_3 = [], [], [],[]
    
    # S_0 =================================
    # print("x", x.shape)
    U_r = pad(x)

    U_0_c = rfft(U_r)

    # First low pass filter
    U_1_c = cdgmm(U_0_c, phi['levels'][0])
    # U_1_c = subsample_fourier(U_1_c, k=2 ** J)
    # U_1_c = subsample_fourier(U_1_c, k=2 ** 0)

    # print("U_1_c", U_1_c.shape)
    S_0 = irfft(U_1_c)
    # print("S_0", S_0.shape)

    S_0 = unpad(S_0, 8)


    out_S_0.append({'coef': S_0,
                    'j': (),
                    'n': (),
                    'theta': ()})
    
    # S_1 =================================
    for n1 in range(len(psi)): # 8 / 16 when J>=2
        j1 = psi[n1]['j']
        theta1 = psi[n1]['theta']

        U_1_c = cdgmm(U_0_c, psi[n1]['levels'][0])
        # if j1 > 0:
        #     U_1_c = subsample_fourier(U_1_c, k=2 ** j1)
        U_1_c = ifft(U_1_c)
        U_1_c = modulus(U_1_c)
        U_1_c = rfft(U_1_c)

        # Second low pass filter
        # S_1_c = cdgmm(U_1_c, phi['levels'][j1])
        S_1_c = cdgmm(U_1_c, phi['levels'][0])
        # S_1_c = subsample_fourier(S_1_c, k=2 ** (J - j1))
        S_1_c = subsample_fourier(S_1_c, k=2 ** 1) # first
    
        S_1_r = irfft(S_1_c)
        S_1_r = unpad(S_1_r,4)

        out_S_1.append({'coef': S_1_r,
                        'j': (j1,),
                        'n': (n1,),
                        'theta': (theta1,)})
        
        # S_2 =================================
        if max_order < 2:
            continue
        for n2 in range(len(psi)):
            j2 = psi[n2]['j']
            theta2 = psi[n2]['theta']

            if j2 <= j1:
                continue
            # only when j2>j1

            # U_2_c = cdgmm(U_1_c, psi[n2]['levels'][j1])
            U_2_c = cdgmm(U_1_c, psi[n2]['levels'][0]) # psi have only one level "0"
            U_2_c = ifft(U_2_c)
            U_2_c = modulus(U_2_c)
            U_2_c = rfft(U_2_c)

            # Third low pass filter
            S_2_c = subsample_fourier(U_2_c, k=2 ** (1))
            S_2_c = cdgmm(S_2_c, phi['levels'][1])
            S_2_c = subsample_fourier(S_2_c, k=2 ** (1))

            S_2_r = irfft(S_2_c)
            S_2_r = unpad(S_2_r,2)

            out_S_2.append({'coef': S_2_r,
                            'j': (j1, j2),
                            'n': (n1, n2),
                            'theta': (theta1, theta2)})
                            
            # S_3 =================================
            for n3 in range(len(psi)):
                j3 = psi[n3]['j']
                theta3 = psi[n3]['theta']

                if j3 <= j2:
                    continue    
                # only when j3>j2>j1

                U_3_c = cdgmm(U_2_c, psi[n3]['levels'][0]) # psi have only one level "0"
                U_3_c = subsample_fourier(U_3_c, k=2 ** 2)
                U_3_c = ifft(U_3_c)
                U_3_c = modulus(U_3_c)
                U_3_c = rfft(U_3_c)

                # Third low pass filter
                S_3_c = cdgmm(U_3_c, phi['levels'][2])
                S_3_c = subsample_fourier(S_3_c, k=2 ** 1)

                S_3_c = irfft(S_3_c)
                S_3_c = unpad(S_3_c,1)

                out_S_3.append({'coef': S_3_c,
                                'j': (j1, j2, j3),
                                'n': (n1, n2, n3),
                                'theta': (theta1, theta2, theta3)})
    # print("len(out_S_0): ", len(out_S_0))
    # print("len(out_S_1): ", len(out_S_1))
    # print("len(out_S_2): ", len(out_S_2))
    # print("len(out_S_3): ", len(out_S_3))


    out_S = []
    out_S.extend(out_S_0)
    out_S.extend(out_S_1)
    out_S.extend(out_S_2)
    out_S.extend(out_S_3)


    if out_type == 'array':
        out_S = concatenate([x['coef'] for x in out_S])

    elif out_type == "ellen_array":
        out_S_0 = concatenate([x['coef'] for x in out_S_0])
        out_S_1 = concatenate([x['coef'] for x in out_S_1])
        out_S_2 = concatenate([x['coef'] for x in out_S_2])
        out_S_3 = concatenate([x['coef'] for x in out_S_3])
        # print("len(out_S_0): ", out_S_0.shape)
        # print("len(out_S_1): ", out_S_1.shape)
        # print("len(out_S_2): ", out_S_2.shape)
        # print("len(out_S_3): ", out_S_3.shape)
        out_S = [out_S_0,out_S_1,out_S_2, out_S_3]

    return out_S


__all__ = ['scattering2d']
