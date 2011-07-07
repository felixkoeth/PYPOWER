# Copyright (C) 1996-2011 Power System Engineering Research Center
# Copyright (C) 2010-2011 Richard Lincoln
#
# PYPOWER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# PYPOWER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PYPOWER. If not, see <http://www.gnu.org/licenses/>.

"""Saves a PYPOWER case file.
"""

from sys import stderr

from numpy import array, c_, r_, any
from scipy.io import savemat

from run_userfcn import run_userfcn

from idx_bus import MU_VMIN, VMIN
from idx_gen import PMIN, MU_PMAX, MU_PMIN, MU_QMIN, MU_QMAX, APF
from idx_brch import MU_ST, MU_SF, BR_STATUS, PF, PT, QT, QF, ANGMAX, MU_ANGMAX
from idx_area import PRICE_REF_BUS
from idx_cost import MODEL, NCOST, PW_LINEAR, POLYNOMIAL


def savecase(fname, ppc, comment=None, version='2'):
    """Saves a PYPOWER case file, given a filename and the data.

    Writes a PYPOWER case file, given a filename and data dict. The C{fname}
    parameter is the name of the file to be created or overwritten. Returns
    the filename, with extension added if necessary. The optional C{comment}
    argument is either string (single line comment) or a list of strings which
    are inserted as comments. When using a PYPOWER case dict, if the
    optional C{version} argument is '1' it will modify the data matrices to
    version 1 format before saving.
    """
    ppc_ver = ppc["version"] = version
    baseMVA, bus, gen, branch = \
        ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]
    areas = ppc["areas"] if "areas" in ppc else None
    gencost = ppc["gencost"] if "gencost" in ppc else None

    ## modifications for version 1 format
    if ppc_ver == "1":
        raise NotImplementedError
#        ## remove extra columns of gen
#        if gen.shape[1] >= MU_QMIN:
#            gen = c_[gen[:, :PMIN], gen[:, MU_PMAX:MU_QMIN]]
#        else:
#            gen = gen[:, :PMIN]
#        ## use the version 1 values for column names
#        shift = MU_PMAX - PMIN - 1
#        tmp = array([MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN]) - shift
#        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN = tmp
#
#        ## remove extra columns of branch
#        if branch.shape[1] >= MU_ST:
#            branch = c_[branch[:, :BR_STATUS], branch[:, PF:MU_ST]]
#        elif branch.shape[1] >= QT:
#            branch = c_[branch[:, :BR_STATUS], branch[:, PF:QT]]
#        else:
#            branch = branch[:, :BR_STATUS]
#        ## use the version 1 values for column names
#        shift = PF - BR_STATUS - 1
#        tmp = array([PF, QF, PT, QT, MU_SF, MU_ST]) - shift
#        PF, QF, PT, QT, MU_SF, MU_ST = tmp

    ## verify valid filename
    l = len(fname)
    rootname = ""
    if l > 2:
        if fname[-3:] == ".py":
            rootname = fname[:-3]
            extension = ".py"
        elif l > 4:
            if fname[-4:] == ".mat":
                rootname = fname[:-4]
                extension = ".mat"

    if not rootname:
        rootname = fname
        extension = ".py"
        fname = rootname + extension

    ## open and write the file
    if extension == ".mat":     ## MAT-file
        savemat(fname, ppc)
    else:                       ## Python file
        try:
            fd = open(fname, "wb")
        except Exception, detail:
            stderr.write("savecase: %s.\n" % detail)
            return fname

        ## function header, etc.
        if ppc_ver == "1":
            if (areas != None) and (gencost != None) and (len(gencost) > 0):
                fd.write('function [baseMVA, bus, gen, branch, areas, gencost] = %s\n' % rootname)
            else:
                fd.write('function [baseMVA, bus, gen, branch] = %s\n' % rootname)
            prefix = ''
        else:
            fd.write('def %s():\n' % rootname)
            prefix = 'ppc'
        if comment:
            if isinstance(comment, basestring):
                fd.write('#%s\n' % comment)
            elif isinstance(comment, list):
                for c in comment:
                    fd.write('#%s\n' % c)
        fd.write('\n#### PYPOWER Case Format : Version %s\n' % ppc_ver)
        if ppc_ver != "1":
            fd.write('ppc[\'version\'] = \'%s\'\n' % ppc_ver)
        fd.write('\n####-----  Power Flow Data  -----####\n')
        fd.write('#### system MVA base\n')
        fd.write('%s[\'sbaseMVA\'] = %g\n' % (prefix, baseMVA))

        ## bus data
        ncols = bus.shape[1]
        fd.write('\n#### bus data\n')
        fd.write('## bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin')
        if ncols >= MU_VMIN + 1:             ## opf SOLVED, save with lambda's & mu's
            fd.write('lam_P lam_Q mu_Vmax mu_Vmin')
        fd.write('\n%s[\'bus\'] = array([\n' % prefix)
        if ncols < MU_VMIN + 1:              ## opf NOT SOLVED, save without lambda's & mu's
            for i in range(bus.shape[0]):
                fd.write('[%d, %d, %g, %g, %g, %g, %d, %.8g, %.8g, %g, %d, %g, %g],\n' % tuple(bus[i, :VMIN + 1]))
        else:                            ## opf SOLVED, save with lambda's & mu's
            for i in range(bus.shape[0]):
                fd.write('[%d, %d, %g, %g, %g, %g, %d, %.8g, %.8g, %g, %d, %g, %g, %.4f, %.4f, %.4f, %.4f],\n' % tuple(bus[:, :MU_VMIN + 1]))
        fd.write('])\n')

        ## generator data
        ncols = gen.shape[1]
        fd.write('\n#### generator data\n')
        fd.write('## bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin')
        if ppc_ver != "1":
            fd.write(' Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q apf')
        if ncols >= MU_QMIN + 1:             # opf SOLVED, save with mu's
            fd.write(' mu_Pmax mu_Pmin mu_Qmax mu_Qmin')
        fd.write('\n%s[\'gen\'] = array([\n' % prefix)
        if ncols < MU_QMIN + 1:              ## opf NOT SOLVED, save without mu's
            if ppc_ver == "1":
                for i in range(gen.shape[0]):
                    fd.write('[%d, %g, %g, %g, %g, %.8g, %g, %d, %g, %g],\n' % tuple(gen[i, :PMIN + 1]))
            else:
                for i in range(gen.shape[0]):
                    fd.write('[%d, %g, %g, %g, %g, %.8g, %g, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g],\n' % tuple(gen[i, :APF + 1]))
        else:
            if ppc_ver == "1":
                for i in range(gen.shape[0]):
                    fd.write('[%d, %g, %g, %g, %g, %.8g, %g, %d, %g, %g, %.4f, %.4f, %.4f, %.4f],\n' % tuple(gen[i, :MU_QMIN + 1]))
            else:
                for i in range(gen.shape[0]):
                    fd.write('[%d, %g, %g, %g, %g, %.8g, %g, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %.4f, %.4f, %.4f, %.4f],\n' % tuple(gen[i, :MU_QMIN + 1]))
        fd.write('])\n')

        ## branch data
        ncols = branch.shape[1]
        fd.write('\n#### branch data\n')
        fd.write('## fbus tbus r x b rateA rateB rateC ratio angle status')
        if ppc_ver != "1":
            fd.write(' angmin angmax')
        if ncols >= QT + 1:                  ## power flow SOLVED, save with line flows
            fd.write(' Pf Qf Pt Qt')
        if ncols >= MU_ST + 1:               ## opf SOLVED, save with mu's
            fd.write(' mu_Sf mu_St')
            if ppc_ver != "1":
                fd.write(' mu_angmin mu_angmax')
        fd.write('\n%s[\'branch\'] = array([\n' % prefix)
        if ncols < QT + 1:                   ## power flow NOT SOLVED, save without line flows or mu's
            if ppc_ver == "1":
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d],\n' % tuple(branch[i, :BR_STATUS + 1]))
            else:
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d, %g, %g],\n' % tuple(branch[i, :ANGMAX + 1]))
        elif ncols < MU_ST + 1:            ## power flow SOLVED, save with line flows but without mu's
            if ppc_ver == "1":
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d, %.4f, %.4f, %.4f, %.4f],\n' % tuple(branch[i, :QT + 1]))
            else:
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d, %g, %g, %.4f, %.4f, %.4f, %.4f],\n' % tuple(branch[i, :QT + 1]))
        else:                            ## opf SOLVED, save with lineflows & mu's
            if ppc_ver == "1":
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f],\n' % tuple(branch[i, :MU_ST + 1]))
            else:
                for i in range(branch.shape[0]):
                    fd.write('[%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %d, %g, %g, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f],\n' % tuple(branch[i, :MU_ANGMAX + 1]))
        fd.write('];\n')

        ## OPF data
        if (areas != None) and (len(areas) > 0) or (gencost != None) and (len(gencost) > 0):
            fd.write('\n####-----  OPF Data  -----####')
        if (areas != None) and (len(areas) > 0):
            ## area data
            fd.write('\n#### area data\n')
            fd.write('## area refbus\n')
            fd.write('%s[\'areas\'] = array([\n' % prefix)
            if len(areas) > 0:
                for i in range(areas.shape[0]):
                    fd.write('[%d, %d],\n' % tuple(areas[i, :PRICE_REF_BUS + 1]))
            fd.write('])\n')
        if gencost != None and len(gencost) > 0:
            ## generator cost data
            fd.write('\n#### generator cost data\n')
            fd.write('## 1 startup shutdown n x1 y1 ... xn yn\n')
            fd.write('## 2 startup shutdown n c(n-1) ... c0\n')
            fd.write('%s[\'gencost\'] = array([\n' % prefix)
            if len(gencost > 0):
                if any(gencost[:, MODEL] == PW_LINEAR):
                    n1 = 2 * max(gencost[gencost[:, MODEL] == PW_LINEAR,  NCOST])
                else:
                    n1 = 0
                if any(gencost[:, MODEL] == POLYNOMIAL):
                    n2 =     max(gencost[gencost[:, MODEL] == POLYNOMIAL, NCOST])
                else:
                    n2 = 0
                n = int( max([n1, n2]) )
                if gencost.shape[1] < n + 4:
                    stderr.write('savecase: gencost data claims it has more columns than it does\n')
                template = '[%d, %g, %g, %d],'
                for i in range(n):
                    template = template + ', %g'
                template = template + '],\n'
                for i in range(gencost.shape[0]):
                    fd.write(template % tuple(gencost[i]))
            fd.write('])\n')

        ## generalized OPF user data
        if ("A" in ppc) and (len(ppc["A"]) > 0) or ("N" in ppc) and (len(ppc["N"]) > 0):
            fd.write('\n####-----  Generalized OPF User Data  -----####')

        ## user constraints
        if ("A" in ppc) and (len(ppc["A"]) > 0):
            ## A
            fd.write('\n#### user constraints\n')
            print_sparse(fd, prefix + "['A']", ppc["A"])
            if ("l" in ppc) and (len(ppc["l"]) > 0) and ("u" in ppc) and (len(ppc["u"]) > 0):
                fd.write('lu = array([\n')
                for i in range(len(l)):
                    fd.write('[%g, %g],\n' % (ppc["l"][i], ppc["u"][i]))
                fd.write('])\n')
                fd.write('%s[\'l\'] = lu[:, 0]\n' % prefix)
                fd.write('%s[\'u\'] = lu[:, 1]\n\n', prefix)
            elif ("l" in ppc) and (len(ppc["l"]) > 0):
                fd.write('%s[\'l\'] = array([\n' % prefix)
                for i in range(len(l)):
                    fd.write('[%g],\n', ppc["l"][i])
                fd.write('])\n\n')
            elif ("u" in ppc) and (len(ppc["u"]) > 0):
                fd.write('%s[\'u\'] = array([\n' % prefix)
                for i in range(len(l)):
                    fd.write('[%g],\n', ppc["u"][i])
                fd.write('])\n\n')

        ## user costs
        if ("N" in ppc) and (len(ppc["N"]) > 0):
            fd.write('\n#### user costs\n')
            print_sparse(fd, prefix + "['N']", ppc["N"])
            if ("H" in ppc) and (len(ppc["H"]) > 0):
                print_sparse(fd, prefix + "['H']", ppc["H"])
            if ("fparm" in ppc) and (len(ppc["fparm"]) > 0):
                fd.write('Cw_fparm = array([\n')
                for i in range(ppc["Cw"]):
                    fd.write('[%g, %d, %g, %g, %g],\n' % tuple(ppc["Cw"][i]) + tuple(ppc["fparm"][i, :]))
                fd.write('])\n')
                fd.write('%s[\'Cw\']    = Cw_fparm[:, 0]\n' % prefix)
                fd.write('%s[\'fparm\'] = Cw_fparm[:, 1:5]\n' % prefix)
            else:
                fd.write('%s[\'Cw\'] = array([\n', prefix)
                for i in range(len(ppc["Cw"])):
                    fd.write('[%g],\n' % ppc["Cw"][i])
                fd.write('])\n')

        ## user vars
        if ('z0' in ppc) or ('zl' in ppc) or ('zu' in ppc):
            fd.write('\n#### user vars\n')
            if ('z0' in ppc) and (len(ppc['z0']) > 0):
                fd.write('%["z0"] = array([\n' % prefix)
                for i in range(len(ppc['z0'])):
                    fd.write('[%g],\n' % ppc["z0"])
                fd.write('])\n')
            if ('zl' in ppc) and (len(ppc['zl']) > 0):
                fd.write('%s["zl"] = array([\n' % prefix)
                for i in range(len(ppc['zl'])):
                    fd.write('[%g],\n' % ppc["zl"])
                fd.write('])\n')
            if ('zu' in ppc) and (len(ppc['zu']) > 0):
                fd.write('%s["zu"] = array([\n' % prefix)
                for i in range(len(ppc['zu'])):
                    fd.write('[%g],\n' % ppc["zu"])
                fd.write('])\n')

        ## execute userfcn callbacks for 'savecase' stage
        if 'userfcn' in ppc:
            run_userfcn(ppc["userfcn"], 'savecase', ppc, fd, prefix)

        ## close file
        fd.close()

    return fname


def print_sparse(fd, varname, A):
    A = A.tocoo()
    i, j, s = A.row, A.col, A.data
    m, n = A.shape

    if len(s) == 0:
        fd.write('%s = sparse((%d, %d))\n' % (varname, m, n))
    else:
        fd.write('ijs = array([\n')
    for k in range(len(i)):
        fd.write('[%d, %d, %g],\n' % (i[k], j[k], s[k]))

    fd.write('])\n')
    fd.write('%s = sparse(ijs[:, 0], ijs[:, 1], ijs[:, 2], %d, %d)\n' % (varname, m, n))
