# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Runs a power flow.
"""

from sys import stdout, stderr

from os.path import dirname, join

from time import time

from numpy import r_, c_, ix_, zeros, pi, ones, exp, argmax,angle
from numpy import flatnonzero as find

#from pypower.bustypes import bustypes
#from pypower.ext2int import ext2int
#from pypower.loadcase import loadcase
#from pypower.ppoption import ppoption
#from pypower.ppver import ppver
#from pypower.makeBdc import makeBdc
from pypower.makeSbus import makeSbus
#from pypower.dcpf import dcpf
#from pypower.makeYbus import makeYbus
from pypower.newtonpf_fast import newtonpf_fast
#from pypower.fdpf import fdpf
#from pypower.gausspf import gausspf
#from pypower.makeB import makeB
#from pypower.pfsoln import pfsoln
#from pypower.printpf import printpf
#from pypower.savecase import savecase
#from pypower.int2ext import int2ext

from pypower.idx_bus import PD, QD, VM, VA, GS, BUS_TYPE, PQ, REF
from pypower.idx_brch import PF, PT, QF, QT
from pypower.idx_gen import PG, QG, VG, QMAX, QMIN, GEN_BUS, GEN_STATUS


def runpf_fast(Ybus, Yf,Yt,ref, pv, pq,on,ppc, ppopt=None, fname='', solvedcase=''):
    """Runs a power flow.

    Runs a power flow [full AC Newton's method by default] and optionally
    returns the solved values in the data matrices, a flag which is C{True} if
    the algorithm was successful in finding a solution, and the elapsed
    time in seconds. All input arguments are optional. If C{casename} is
    provided it specifies the name of the input data file or dict
    containing the power flow data. The default value is 'case9'.

    If the ppopt is provided it overrides the default PYPOWER options
    vector and can be used to specify the solution algorithm and output
    options among other things. If the 3rd argument is given the pretty
    printed output will be appended to the file whose name is given in
    C{fname}. If C{solvedcase} is specified the solved case will be written
    to a case file in PYPOWER format with the specified name. If C{solvedcase}
    ends with '.mat' it saves the case as a MAT-file otherwise it saves it
    as a Python-file.

    If the C{ENFORCE_Q_LIMS} options is set to C{True} [default is false] then
    if any generator reactive power limit is violated after running the AC
    power flow, the corresponding bus is converted to a PQ bus, with Qg at
    the limit, and the case is re-run. The voltage magnitude at the bus
    will deviate from the specified value in order to satisfy the reactive
    power limit. If the reference bus is converted to PQ, the first
    remaining PV bus will be used as the slack bus for the next iteration.
    This may result in the real power output at this generator being
    slightly off from the specified values.

    Enforcing of generator Q limits inspired by contributions from Mu Lin,
    Lincoln University, New Zealand (1/14/05).

    @author: Ray Zimmerman (PSERC Cornell)
    """
    ## default arguments


    ## options



    ## read data
    #ppc = loadcase(casedata)
    
    ## convert to internal indexing
    

    ppc["branch"][:,[0,1]]-=1
    ppc["bus"][:,0]-=1
    ppc["gen"][:,0]-=1

    baseMVA, bus, gen, branch = \
        ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]

    ## get bus index lists of each type of bus
    #ref, pv, pq = bustypes(bus, gen)

    #
    # generator info
    #print(gen[:, GEN_STATUS])
    #on = find(gen[:, GEN_STATUS] > 0)      ## which generators are on?
    gbus = gen[on, GEN_BUS].astype(int)    ## what buses are they at?

    ##-----  run the power flow  -----
    t0 = time()


    V0  = bus[:, VM] * exp(1j * 0.017453292519943295 * bus[:, VA])
    V0[gbus] = gen[on, VG] / abs(V0[gbus]) * V0[gbus]


    ## build admittance matrices
    #Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)

    ## compute complex bus power injections [generation - load]
    Sbus = makeSbus(baseMVA, bus, gen)

    ## run the power flow

    V, success, i = newtonpf_fast(Ybus, Sbus, V0, ref, pv, pq, ppopt)

    ## update data matrices with solution
    #bus, gen, branch = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq)
    bus[:, VM] = abs(V)
    bus[:, VA] = angle(V) * 180 / pi

    #UNTIL HERE
    ppc["et"] = time() - t0
    ppc["success"] = success

    ##-----  output results  -----
    ## convert back to original bus numbering & print results
    ppc["bus"], ppc["gen"], ppc["branch"] = bus, gen, branch

    ppc["branch"][:,[0,1]]+=1
    ppc["bus"][:,0]+=1
    ppc["gen"][:,0]+=1

    return ppc, success,i

if __name__ == '__main__':
    runpf()
