import mysql.connector
from pyscf import gto, dft,scf,cc
import  math
import numpy as np
import numpy
from pyscf.geomopt.geometric_solver import (optimize)

def lyp(a,b,gaa,gbb,gnn,i):
    A = 0.06513211737123356
    B = 0.16550543387395644
    C = 0.1694065214864509
    Dd = 0.5807106359304246
    n=a+b
    CF = (3.0 / 10.0) * math.pow(3.0 * math.pow(math.pi, 2), 2.0 / 3.0)

    icbrtn = n** (-1.0 / 3.0)
    P = 1 / (1 + Dd * icbrtn)
    omega = numpy.exp(-C * icbrtn) * P * n** (-11.0 / 3.0)
    delta = icbrtn * (C + Dd * P)
    n2 = n * n
    f = -A * (4 * a * b * P / n
              + B * omega
              * (a * b
                 * (math.pow(2, 11.0 / 3.0) * CF
                    * ((a** (8.0 / 3.0)) + (b**( 8.0 / 3.0)))
                    + (47.0 - 7.0 * delta) * gnn / 18.0
                    - (2.5 - delta / 18.0) * (gaa + gbb)
                    - (delta - 11.0) / 9.0 * (a * gaa + b * gbb) / n)
                 - 2.0 / 3.0 * n2 * gnn + (2.0 / 3.0 * n2 - a * a) * gbb
                 + (2.0 / 3.0 * n2 - b * b) * gaa))
    return f



def eval_xc_lyp_b(xc_code, rho, spin, relativity=0, deriv=2, verbose=None,omega=None):
    rho1 = rho[0]
    rho2 = rho[1]

    a, dx1, dy1, dz1 = rho1[:4]
    b, dx2, dy2, dz2 = rho2[:4]
    #a=a+1E-20
    #b = b + 1E-20

    gaa = dx1 ** 2 + dy1 ** 2 + dz1 ** 2
    gbb = dx2 ** 2 + dy2 ** 2 + dz2 ** 2
    gnn = (dx1+dx2) ** 2 + (dy1+dy2) ** 2 + (dz1+dz2) ** 2

    delta1=a*1.0E-9
    delta2 =b*1.0E-9

    exc=lyp(a,b,gaa,gbb,gnn,1)
    exc1 = (lyp(a+delta1, b, gaa, gbb, gnn,2))
    exc2 = lyp(a , b+delta2, gaa, gbb, gnn,3)

    vrho1=(exc1 - exc) / (delta1+1e-200)
    vrho2 = (exc2 - exc) / (delta2+1e-200)

    delta1 = (gaa)*1.0E-25
    delta2 = (gbb ) * 1.0E-25
    delta3 = (gnn ) * 1.0E-25


    exc1g = lyp(a , b, gaa+delta1, gbb, gnn,4)
    exc2g = lyp(a, b , gaa, gbb+delta2, gnn,5)
    exc3g = lyp(a, b, gaa, gbb , gnn+delta3,6)

    vgama1 = (exc1g - exc) / (delta1+1e-250)
    vgama2 = (exc2g - exc) / (delta2+1e-250)
    vgama3 = (exc3g - exc) / (delta3+1e-250)

    exc1=exc1/(a+b)
    exc2 = exc2 / (a + b)



    b88 = dft.libxc.eval_xc('B88,', rho, spin, relativity, deriv,
                               verbose)


    vgamma_ = np.array([vgama1, vgama3, vgama2])
    vgamma_ = vgamma_ + np.array(b88[1][1]).T

    vgamma = np.transpose(vgamma_)



    vrho_ = np.array([vrho1, vrho2])
    vrho_=vrho_+numpy.array(b88[1][0]).T


    vrho = np.transpose(vrho_)

    exc1 = np.array([exc1])
    exc2 = np.array([exc2])

    exc = np.transpose(0.5*exc1 + 0.5*exc2+b88[0])
    vxc = (vrho, vgamma, None, None)
    fxc = None  # 2nd order functional derivative
    kxc = None  # 3rd order functional derivative
    return exc, vxc, fxc, kxc



mol = gto.Mole()
mol.verbose = 4
mol.atom = """8	0.0000	0.0000	0.1173
1	0.0000	0.7572	-0.4692
1	0.0000	-0.7572	-0.4692"""
mol.charge = 0
mol.spin = 0
mol.basis = "6-31G"
mol.build()
mf = scf.UHF(mol)
mol = optimize(mf, maxsteps=100)
mfl = dft.UKS(mol)
mfl = mfl.define_xc_(eval_xc_lyp_b, xctype='GGA')
mfl.level_shift = 0.5
mfl.grids.level = 0
mfl.conv_tol = 0.000000001
ae =  mfl.kernel()
print(ae)