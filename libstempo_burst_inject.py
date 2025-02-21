import numpy as N
import ephem

def add_burst(psr, gwtheta, gwphi, waveform_plus, waveform_cross, psi=0.0, tref=0, remove_quad=False):
    """
    Function to create GW-induced residuals from an arbitrary GW waveform assuming elliptical polarization.
    :param psr: pulsar object
    :param gwtheta: Polar angle of GW source in celestial coords [radians]
    :param gwphi: Azimuthal angle of GW source in celestial coords [radians]
    :param waveform_plus: Function defining the plus polarized waveform of the GW [function]
    :param waveform_cross: Function defining the cross polarized waveform of the GW [function]
    :param psi: Polarization angle [radians]. Mixes h+ and hx corresponding to rotation along the propagation direction (see eq. (7.24-25) in Maggiore Vol1, 2008).
    :param tref: Start time, such that gw_waveform gets t-tref as the time argument
    :param remove_quad: Fit out quadratic from residual if True to simulate f and fdot timing fit.
    :returns: Vector of induced residuals
    """

    # define variable for later use
    cosgwtheta, cosgwphi = N.cos(gwtheta), N.cos(gwphi)
    singwtheta, singwphi = N.sin(gwtheta), N.sin(gwphi)

    # unit vectors to GW source
    m = N.array([singwphi, -cosgwphi, 0.0])
    n = N.array([-cosgwtheta*cosgwphi, -cosgwtheta*singwphi, singwtheta])
    omhat = N.array([-singwtheta*cosgwphi, -singwtheta*singwphi, -cosgwtheta])

    # pulsar location
    if 'RAJ' and 'DECJ' in psr.pars():
        ptheta = N.pi/2 - psr['DECJ'].val
        pphi = psr['RAJ'].val
    elif 'ELONG' and 'ELAT' in psr.pars():
        fac = 180./N.pi
        coords = ephem.Equatorial(ephem.Ecliptic(str(psr['ELONG'].val*fac),
                                                 str(psr['ELAT'].val*fac)))

        ptheta = N.pi/2 - float(repr(coords.dec))
        pphi = float(repr(coords.ra))

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = N.array([N.sin(ptheta)*N.cos(pphi), N.sin(ptheta)*N.sin(pphi),\
            N.cos(ptheta)])

    #print(ptheta, pphi, gwtheta, gwphi)

    fplus = 0.5 * (N.dot(m, phat)**2 - N.dot(n, phat)**2) / (1+N.dot(omhat, phat))
    fcross = (N.dot(m, phat)*N.dot(n, phat)) / (1 + N.dot(omhat, phat))

    # get toas from pulsar object
    toas = psr.toas()*86400 - tref

    # define residuals: hplus and hcross
    hplus = waveform_plus(toas)
    hcross = waveform_cross(toas)

    #apply rotation by psi angle (see e.g. eq. (7.24-25) in Maggiore Vol1, 2008)
    rplus = hplus*N.cos(2*psi) - hcross*N.sin(2*psi)
    rcross = hplus*N.sin(2*psi) + hcross*N.cos(2*psi)

    # residuals
    res = -fplus*rplus - fcross*rcross

    if remove_quad:
        pp = N.polyfit(N.array(toas, dtype=N.double), N.array(res, dtype=N.double), 2)
        res = res - pp[0]*toas**2 -pp[1]*toas - pp[2]

    psr.stoas[:] += res/86400

    return res

def add_glitch(psr, waveform, tref=0):
    """
    Function to create incoherent residuals of arbitrary waveform in a given pulsar.
    :param psr: pulsar object
    :param waveform: Function defining the waveform of the glitch [function]
    :param tref: Start time, such that gw_waveform gets t-tref as the time argument
    :returns: Vector of induced residuals
    """

    # get toas from pulsar object
    toas = psr.toas()*86400 - tref

    # call waveform function
    res = waveform(toas)

    psr.stoas[:] += res/86400
