#imports
import os, glob, json
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sps

from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import selections
from enterprise.signals.selections import Selection

from enterprise_extensions import models, hypermodel
from enterprise_extensions.model_utils import bayes_fac
from la_forge import core as co
from la_forge import diagnostics as dg

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

import pickle
import time

pkl_file = '/home/mitch/test_quick_burst/15_year_data/psrs_trimmed_SNR99p.pkl'
with open(pkl_file, 'rb') as file:
    psrs = pickle.load(file)

parent_dir = '/home/mitch/test_quick_burst/15_year_data/15_year_CURN_run/delete/'

psr = psrs[17]
selection = selections.Selection(selections.by_backend)
selection_ng = selections.Selection(selections.nanograv_backends)

#WN
efac = parameter.Uniform(0.01,10.0)
equad = parameter.Uniform(-8.5,-5)
ecorr = parameter.Uniform(-8.5,-5)

wn = white_signals.MeasurementNoise(efac=efac, log10_t2equad=equad)
ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection_ng)

#timing model
tm = gp_signals.MarginalizingTimingModel(use_svd=True)

#RN
log10_A = parameter.Uniform(-20,-11)("gw_curn_amp")
gamma = parameter.Uniform(0,7)('gw_gamma')

tmax = psr.toas.max()
tmin = psr.toas.min()
tspan = (tmax-tmin)
bins = 30
freq = np.arange(1/tspan, bins/tspan, 1/tspan)
pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
rn = gp_signals.FourierBasisGP(spectrum=pl, modes=freq)

#full model
s = tm+ wn + ec + rn
pta = signal_base.PTA(s(psr))
x0 = np.hstack([p.sample() for p in pta.params])
ndim = len(x0)

cov = np.diag(np.ones(ndim)*0.01**2)

output_dir = parent_dir + psr.name
os.makedirs(output_dir, exist_ok=True)

sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, outDir=output_dir, resume=False)
print(pta.params)
sampler_steps = 10e5

time_start = time.perf_counter()
sampler.sample(x0,sampler_steps)
time_end = time.perf_counter()
print("elapsed time: {}".format(time_end-time_start))




