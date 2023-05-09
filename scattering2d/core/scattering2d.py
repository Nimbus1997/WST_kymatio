# modify 2
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
    out_S_0, out_S_1, out_S_2 = [], [], []
    
    print("x", x.shape)
    U_r = pad(x)
    print("U_r", U_r.shape)

    U_0_c = rfft(U_r)
    print("U_0_c", U_0_c.shape)



    # First low pass filter
    U_1_c = cdgmm(U_0_c, phi['levels'][0])
    print("U_1_c-cdgmm", U_1_c.shape)

    # U_1_c = subsample_fourier(U_1_c, k=2 ** J) # original
    U_1_c = subsample_fourier(U_1_c, k=2 ** 0)
    print("U_1_c-subsampelf", U_1_c.shape)


    S_0 = irfft(U_1_c)
    print("S_0", S_0.shape)

    S_0 = unpad(S_0)
    print("S_0", S_0.shape)
    print("--------1----------")


    out_S_0.append({'coef': S_0,
                    'j': (),
                    'n': (),
                    'theta': ()})

    for n1 in range(len(psi)): # 8 / 16 when J>=2
        j1 = psi[n1]['j']
        theta1 = psi[n1]['theta']
        print("U_0_c-start", U_0_c.shape)
        U_1_c = cdgmm(U_0_c, psi[n1]['levels'][0])
        print("U_1_c-cdgmm", U_1_c.shape)

        if j1 > 0:
            # U_1_c = subsample_fourier(U_1_c, k=2 ** j1) # original
            U_1_c = subsample_fourier(U_1_c, k
                                      =2 ** 1)
            print("U_1_c-subsample", U_1_c.shape)

        U_1_c = ifft(U_1_c)
        print("U_1_c-ifft", U_1_c.shape)

        U_1_c = modulus(U_1_c)
        print("U_1_c-modulus", U_1_c.shape)

        U_1_c = rfft(U_1_c)
        print("U_1_c-rfft", U_1_c.shape)


        # Second low pass filter
        # S_1_c = cdgmm(U_1_c, phi['levels'][j1])
        # print("S_1_c-cdgmm", S_1_c.shape)

        # S_1_c = subsample_fourier(S_1_c, k=2 ** (J - j1)) #original
        S_1_c = subsample_fourier(U_1_c, k=2 ** 1)
        print("S_1_c-subsample", S_1_c.shape)


        S_1_r = irfft(S_1_c)
        print("S_1_r-irfft", S_1_r.shape)

        S_1_r = unpad(S_1_r)
        print("S_1_r-unpad", S_1_r.shape)
        
        out_S_1.append({'coef': S_1_r,
                        'j': (j1,),
                        'n': (n1,),
                        'theta': (theta1,)})

        if max_order < 2:
            continue
        print("-----------2-----------")
        for n2 in range(len(psi)):
            j2 = psi[n2]['j']
            theta2 = psi[n2]['theta']
            if j2 <= j1:
                continue
            print("U_1_c-start", U_1_c.shape)
            U_2_c = cdgmm(U_1_c, psi[n2]['levels'][j1])
            print("U_2_c-cdgmm", U_2_c.shape)

            # U_2_c = subsample_fourier(U_2_c, k=2 ** (j2 - j1)) # original
            U_2_c = subsample_fourier(U_2_c, k=2 ** 1)
            print("U_2_c-subsample", U_2_c.shape)

            U_2_c = ifft(U_2_c)
            print("U_2_c-ifft", U_2_c.shape)

            U_2_c = modulus(U_2_c)
            print("U_2_c-modulus", U_2_c.shape)
            U_2_c = rfft(U_2_c)
            print("U_2_c-rfft", U_2_c.shape)


            # Third low pass filter
            # S_2_c = cdgmm(U_2_c, phi['levels'][j2]) # 여기서 문제일어남 ! ()
            # print("S_2_c-cdgmm", S_2_c.shape)

            # S_2_c = subsample_fourier(S_2_c, k=2 ** (J - j2)) #original
            S_2_c = subsample_fourier(U_2_c, k=2 **1)
            print("S_2_c-subsample", S_2_c.shape)

            S_2_r = irfft(S_2_c)
            print("S_2_r-irfft", S_2_r.shape)

            S_2_r = unpad(S_2_r)
            print("S_2_r-unpad", S_2_r.shape)


            out_S_2.append({'coef': S_2_r,
                            'j': (j1, j2),
                            'n': (n1, n2),
                            'theta': (theta1, theta2)})

    out_S = []

    out_S.extend(out_S_0)
    out_S.extend(out_S_1)
    out_S.extend(out_S_2)

    if out_type == 'array':
        out_S = concatenate([x['coef'] for x in out_S])

    return out_S


__all__ = ['scattering2d']
