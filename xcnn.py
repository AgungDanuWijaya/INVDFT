from pyscf import gto, dft,scf,cc
import  math
import numpy as np
import numpy
from pyscf.geomopt.geometric_solver import (optimize)



def sigmoid(x):
    return x

def calcsig(h):
    return [sigmoid(x) for x in h]

def calc(input, w, p):
    hasil = [0] * len(w)
    for i in range(len(w)):
        hasil[i]=w[i] * np.power(np.abs(input), p[i])
        #print(hasil[i],w[i],input,p[i])

    return hasil

def f_ANN(input, c_node, w):
    node = c_node[1:]
    for k in range(len(node)):
        re = [0] * node[k]
        for i in range(len(input)):
            re = numpy.array(re)+numpy.array( calc(input[i], w[0][k][i], w[1][k][i]))
            #print(re)
        input = calcsig(re)
    return np.sum(np.array(input))


def annn(rho,gama) :
    c_node = [2, 4, 4, 1]

    we = [0] * (len(c_node) - 1)
    for i in range(len(c_node) - 1):
        re_i = [0] * c_node[i]
        for j in range(c_node[i]):
            re_j = [0] * c_node[i + 1]
            re_i[j] = re_j
        we[i] = re_i
    we1 = [0] * (len(c_node) - 1)
    for i in range(len(c_node) - 1):
        re_i1 = [0] * c_node[i]
        for j in range(c_node[i]):
            re_j1 = [0] * c_node[i + 1]
            re_i1[j] = re_j1
        we1[i] = re_i1
    param = [0, 0]
    w = [we, we1, param]

    w[0][0][0][0] = -0.4642460590404523
    w[0][0][0][1] = 0.03932529165484453
    w[0][0][0][2] = -0.05742808728581002
    w[0][0][0][3] = 0.1448278459670839
    w[0][0][1][0] = -0.8414675419844087
    w[0][0][1][1] = 0.09237642403600035
    w[0][0][1][2] = -0.5350224699956617
    w[0][0][1][3] = 0.16858223408404713
    w[0][1][0][0] = -0.6322469103869556
    w[0][1][0][1] = -0.17100016444625876
    w[0][1][0][2] = 1.4687374018993005
    w[0][1][0][3] = -0.031734908407939105
    w[0][1][1][0] = -0.23091565868261116
    w[0][1][1][1] = -0.9151282211144409
    w[0][1][1][2] = 0.47328275018280386
    w[0][1][1][3] = -0.39008791629439027
    w[0][1][2][0] = -0.5407581178784329
    w[0][1][2][1] = -0.8477526122275791
    w[0][1][2][2] = -0.274803138048219
    w[0][1][2][3] = -0.570666538609607
    w[0][1][3][0] = 0.48808934263130377
    w[0][1][3][1] = 0.9010902662630734
    w[0][1][3][2] = 0.8970276926153301
    w[0][1][3][3] = 0.1590492497274638
    w[0][2][0][0] = -0.13987582098243126
    w[0][2][1][0] = 0.1027814352341956
    w[0][2][2][0] = -0.3827148234703167
    w[0][2][3][0] = 0.08013101914988155
    w[1][0][0][0] = 1.288053810303302
    w[1][0][0][1] = 1.2928371534171186
    w[1][0][0][2] = 1.4601921475108928
    w[1][0][0][3] = 1.5754409008238173
    w[1][0][1][0] = 1.1548457227980655
    w[1][0][1][1] = 0.971073408124568
    w[1][0][1][2] = 1.0982615302008154
    w[1][0][1][3] = 1.1013800243365959
    w[1][1][0][0] = 0.8667318937617526
    w[1][1][0][1] = 0.6303016698137959
    w[1][1][0][2] = 1.1109660903823309
    w[1][1][0][3] = 0.7225172336323346
    w[1][1][1][0] = 1.6018171620385355
    w[1][1][1][1] = 1.0079246160214066
    w[1][1][1][2] = 1.9407476116901141
    w[1][1][1][3] = 0.8096917634625275
    w[1][1][2][0] = 0.8135232363169674
    w[1][1][2][1] = 1.110023178749941
    w[1][1][2][2] = 1.0790625180171292
    w[1][1][2][3] = 1.1735414194276015
    w[1][1][3][0] = 0.7083167642610289
    w[1][1][3][1] = 1.1000294543722615
    w[1][1][3][2] = 1.196790239881579
    w[1][1][3][3] = 1.4268609385411568
    w[1][2][0][0] = 0.9659738022138349
    w[1][2][1][0] = 1.2231412844974938
    w[1][2][2][0] = 0.8524612977316827
    w[1][2][3][0] = 1.0709915478396048
    w[2][0] = 0.9966214012485962
    w[2][1] = 0.10382707907281696


    input = [rho,rho ** w[2][0] * gama ** w[2][1]]
    if (rho < 1E-11):
        input = [0, 0]
    if (gama < 1E-11):
        input = [0, 0]
    if (gama < 1E-11 and rho < 1E-11):
        input = [0, 0]


    a = f_ANN(input, c_node, w)
    return a


def eval_xc_gga(xc_code, rho, spin, relativity=0, deriv=2, verbose=None,omega=None):
    rho1 = rho[0]
    rho2 = rho[1]

    rho01, dx1, dy1, dz1 = rho1[:4]
    rho02, dx2, dy2, dz2 = rho2[:4]
    rho01=rho01+1E-20
    rho02 = rho02 + 1E-20
    w1 = rho01 / (rho01 + rho02)
    w2 = rho02 / (rho01 + rho02)
    gamma1 = dx1 ** 2 + dy1 ** 2 + dz1 ** 2
    gamma2 = dx2 ** 2 + dy2 ** 2 + dz2 ** 2



    exc1=[0]*len(rho01)
    exc2 = [0] * len(rho01)
    vgamma_1=[0] * len(rho01)
    vgamma_3 = [0] * len(rho01)

    vrho1=[0] * len(rho02)
    vrho2 = [0] * len(rho02)
    vgamma_2 = [0] * len(rho02)
    for inm in range(len(rho02)) :
        delta1 = (gamma1[inm] ) * 1.0E-25
        delta2 = (gamma2[inm] ) * 1.0E-25
        delta1_ = (rho01[inm]) * 1.0E-9
        delta2_ = (rho02[inm]) * 1.0E-9

        ex1 = annn(rho01[inm], gamma1[inm])
        ex1_ = annn(rho01[inm] + delta1_, gamma1[inm])
        ex1__ = annn(rho01[inm] , gamma1[inm]+delta1)
        vrho1[inm] = (ex1_ - ex1) / (delta1_+1e-250)
        exc1[inm] = (ex1/rho01[inm])*w1[inm]
        vgamma_1[inm] = (ex1__ - ex1)  / (delta1+1e-250)

        ex2 = annn(rho02[inm], gamma2[inm])
        ex2_ = annn(rho02[inm] + delta2_, gamma2[inm])
        ex2__ = annn(rho02[inm] , gamma2[inm]+delta2)
        vrho2[inm] = (ex2_ - ex2) / (delta2_+1e-250)
        exc2[inm] = (ex2/rho02[inm])*w2[inm]
        vgamma_2[inm] = (ex2__ - ex2) / (delta2+1e-250)

    pbe_xc = dft.libxc.eval_xc(',LYP', rho, spin, relativity, deriv,
                               verbose)
    vgamma_ = np.array([vgamma_1, vgamma_3, vgamma_2])
    vgamma_=vgamma_+np.array(pbe_xc[1][1]).T
    vgamma = np.transpose(vgamma_)
    vrho_ = np.array([vrho1, vrho2])
    vrho_ = vrho_ + numpy.array(pbe_xc[1][0]).T
    vrho = np.transpose(vrho_)
    exc1 = np.array([exc1])
    exc2 = np.array([exc2])
    exc = np.transpose(exc1 + exc2 + pbe_xc[0])


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
mol.basis = "aug-cc-pvdz"
mol.build()
mf = scf.UHF(mol)
mfl = dft.UKS(mol)
mfl = mfl.define_xc_(eval_xc_gga, xctype='GGA')
mfl.level_shift = 0.5
mfl.grids.level = 0
mfl.conv_tol = 0.000000001
ae =  mfl.kernel()
print(ae)
