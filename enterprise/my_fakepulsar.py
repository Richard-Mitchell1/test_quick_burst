import os
import numpy as np
import libstempo

def fakepulsar(parfile, obstimes, toaerr, freq=1440.0, observatory="AXIS", flags="", iters=3):
    """Returns a libstempo tempopulsar object corresponding to a noiseless set
    of observations for the pulsar specified in 'parfile', with observations
    happening at times (MJD) given in the array (or list) 'obstimes', with
    measurement errors given by toaerr (us).
    A new timfile can then be saved with pulsar.savetim(). Re the other parameters:
    - 'toaerr' needs to be either a common error, or a list of errors
       of the same length of 'obstimes';
    - 'freq' can be either a common observation frequency in MHz, or a list;
       it defaults to 1440;
    - 'observatory' can be either a common observatory name, or a list;
       it defaults to the IPTA MDC 'AXIS';
    - 'flags' can be a string (such as '-sys EFF.EBPP.1360') or a list of strings;
       it defaults to an empty string;
    - 'iters' is the number of iterative removals of computed residuals from TOAs
      (which is how the fake pulsar is made...)"""

    import tempfile

    outfile = tempfile.NamedTemporaryFile(delete=False)
    #print(outfile.name)

    outfile.write(b"FORMAT 1\n")
    outfile.write(b"MODE 1\n")

    obsname = "fake_" + os.path.basename(parfile)
    if obsname[-4:] == ".par":
        obsname = obsname[:-4]

    for i, t in enumerate(obstimes):
        tf_line = "{0} {1} {2} {3} {4}".format( obsname, _geti(freq, i), t, _geti(toaerr, i), _geti(observatory, i))
        for key in flags.keys():
            if flags[key][i] != '':
                tf_line += " -{0} {1}".format(key, flags[key][i])
            #else:
            #    tf_line += " -{0} {1}".format(key, "np.nan")
        tf_line += "\n"
        #print(tf_line)
        outfile.write(tf_line.encode("ascii"))

    timfile = outfile.name
    outfile.close()

    pulsar = libstempo.tempopulsar(parfile, timfile, dofit=False, maxobs=60000)

    for i in range(iters):
        pulsar.stoas[:] -= pulsar.residuals() / 86400.0
        pulsar.formbats()

    os.remove(timfile)

    return pulsar


def _geti(x, i):
    return x[i] if isinstance(x, (tuple, list, np.ndarray)) else x
