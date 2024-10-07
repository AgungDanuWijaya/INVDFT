package solid;

import java.io.FileInputStream;
import java.io.FileWriter;   // Import the FileWriter class
import java.io.IOException;  // Import the IOException class to handle errors
import java.io.ObjectInputStream;

public class tulisbobot {
  public static void main(String[] args) throws ClassNotFoundException {
    try {
      FileWriter myWriter = new FileWriter("/home/agung/project/solid/q-e-qe-6.6/Modules/exchange_gga.f90");
    ObjectInputStream inputStream_ = new ObjectInputStream(new FileInputStream("/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/JQC_data/ann_new_dft_monte"));

               double krenl[][][][] = (double[][][][]) inputStream_.readObject();
               String bobot="";
               for (int i = 0; i < krenl.length; i++) {
                   //System.out.println("batas "+i);
                    for (int j = 0; j < krenl[i].length; j++) {
                        for (int k = 0; k < krenl[i][j].length; k++) {
                            for (int l = 0; l < krenl[i][j][k].length; l++) {
                                String kk="";
                                if(i==0){
                                kk="w";
                                }
                                else if(i==1){
                                kk="p";
                                } else if(i==2){
                                kk="n";
                                }
                                //System.out.println("$"+kk+"_{"+j+" "+k+" "+l+" }$="+krenl[i][j][k][l]+""); 
                              bobot+=("w("+(i+1)+","+(j+1)+","+(k+1)+","+(l+1)+")="+krenl[i][j][k][l]+"\n"); 
                            }
                        }

                    }

                }
  bobot+="\n";
      String code1="!\n" +
"MODULE exch_gga !<GPU:exch_gga=>exch_gga_gpu>\n" +
"!\n" +
"CONTAINS\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"\n" +
"function adddot(a, b) result(c)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP), dimension(:), intent(in) :: a, b\n" +
"        REAL(DP), dimension(size(a)) :: c\n" +
"        integer :: i\n" +
"\n" +
"        c = a + b\n" +
"    end function adddot\n" +
"\n" +
"    function sum(a) result(s)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP), dimension(:), intent(in) :: a\n" +
"        REAL(DP) :: s\n" +
"        integer :: i\n" +
"\n" +
"        s = 0.0\n" +
"        do i = 1, size(a)\n" +
"            s = s + a(i)\n" +
"        end do\n" +
"    end function sum\n" +
"\n" +
"    function sigmoid(x) result(s)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP), intent(in) :: x\n" +
"        REAL(DP) :: s\n" +
"\n" +
"        s = x\n" +
"    end function sigmoid\n" +
"\n" +
"    function calcsig(h) result(res)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP), dimension(:), intent(inout) :: h\n" +
"        REAL(DP) :: res(size(h))\n" +
"        integer :: i\n" +
"\n" +
"        do i = 1, size(h)\n" +
"            h(i) = sigmoid(h(i))\n" +
"        end do\n" +
"        res = h\n" +
"    end function calcsig\n" +
"\n" +
"    function calc(input, w, p,ind) result(hasil)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP), intent(in) :: input\n" +
"        integer, intent(in) :: ind\n" +
"        REAL(DP), dimension(:), intent(in) :: w, p\n" +
"        REAL(DP), dimension(size(w)) :: hasil\n" +
"        integer :: i\n" +
"\n" +
"        do i = 1, ind\n" +
"            hasil(i) = w(i) * abs(input)**p(i)\n" +
"        end do\n" +
"    end function calc\n" +
"\n" +
"    function f_ANN(input, c_node, w) result(output)\n" +
"     USE kinds, ONLY: DP\n" +
"        REAL(DP), dimension(:), intent(in) :: input\n" +
"        integer, dimension(:), intent(in) :: c_node\n" +
"        REAL(DP), dimension(:,:,:,:) :: w\n" +
"        integer :: node(size(c_node) - 1)\n" +
"        REAL(DP) :: re(4)\n" +
"        REAL(DP), dimension(4) :: input_\n" +
"        integer :: i, k\n" +
"        REAL(DP)::output\n" +
"\n" +
"        do i = 2, size(c_node)\n" +
"            node(i - 1) = c_node(i)\n" +
"        end do\n" +
"		input_=input\n" +
"        do k = 1, size(node)\n" +
"            re = 0.0\n" +
"            do i = 1, c_node(k)\n" +
"                re = adddot(re, calc(input_(i), w(1, k, i,:), w(2, k, i,:),c_node(k+1)))\n" +
"            end do\n" +
"	    input_=re\n" +
"\n" +
"\n" +
"      \n" +
"        end do\n" +
"        output = input_(1)\n" +
"    end function f_ANN\n" +
"    \n" +
"    function annn(rho, gama)\n" +
"    USE kinds, ONLY: DP\n" +
"        REAL(DP) :: annn\n" +
"        REAL(DP) :: rho, gama\n" +
"        integer, dimension(4) :: c_node\n" +
"        REAL(DP), dimension(3,3,4,4) :: w,w1\n" +
"        REAL(DP), dimension(2) :: input\n" +
"        integer :: i, j, k\n" +
"\n" +
"        c_node = [2, 4, 4, 1]\n" +
"\n" +
"    \n" +
"   "+bobot+
"\n" +
"        \n" +
"        input = [rho, rho**w(3,1,1,1) * gama**w(3,1,1,2)]\n" +
"        if (rho < 1.0E-11) then\n" +
"            input = [0.0, 0.0]\n" +
"        end if\n" +
"        if (gama < 1.0E-11) then\n" +
"            input = [0.0, 0.0]\n" +
"        end if\n" +
"        if (gama < 1.0E-11 .and. rho < 1.0E-11) then\n" +
"            input = [0.0, 0.0]\n" +
"        end if\n" +
"\n" +
"        annn=2.0**(1.0/3.0)*f_ANN(input, c_node,w)\n" +
"\n" +
"    \n" +
"end function annn\n" +
"\n" +
"\n" +
"SUBROUTINE ann_b88( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)\n" +
"  !! only gradient-corrected part, no Slater term included\n" +
"  !\n" +
"  USE kinds, ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee,delta,grho1\n" +
"  REAL(DP), PARAMETER :: beta=0.0042_DP\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, two13=1.259921049894873_DP\n" +
"                                          ! two13= 2^(1/3)\n" +
"  !\n" +
"  \n" +
"grho1=grho\n" +
" delta=rho*1e-9\n" +
"sx=annn(rho,grho1)\n" +
"v1x=(annn(rho+delta,grho1)-sx)/delta\n" +
"delta=grho1*1e-25\n" +
"v2x=(annn(rho,grho1+delta)-sx)/delta\n" +
"\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE ann_b88\n" +
"\n" +
"SUBROUTINE becke88( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"   !-----------------------------------------------------------------------\n" +
"  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)\n" +
"  !! only gradient-corrected part, no Slater term included\n" +
"  !\n" +
"  USE kinds, ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee,delta,grho1\n" +
"  REAL(DP), PARAMETER :: beta=0.0042_DP\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, two13=1.259921049894873_DP\n" +
"                                          ! two13= 2^(1/3)\n" +
"  !\n" +
"\n" +
"\n" +
"grho1=grho*(2**(1.0/3.0))\n" +
" delta=rho*1e-9\n" +
"sx=annn(rho,grho1)\n" +
"v1x=(annn(rho+delta,grho1)-sx)/delta\n" +
"delta=grho1*1e-25\n" +
"v2x=(annn(rho,grho1+delta)-sx)/delta\n" +
"\n" +
" \n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE becke88\n" +
"\n" +
"SUBROUTINE becke88_( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)\n" +
"  !! only gradient-corrected part, no Slater term included\n" +
"  !\n" +
"  USE kinds, ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee\n" +
"  REAL(DP), PARAMETER :: beta=0.0042_DP\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, two13=1.259921049894873_DP\n" +
"                                          ! two13= 2^(1/3)\n" +
"  !\n" +
"  rho13 = rho**third\n" +
"  rho43 = rho13**4\n" +
"  !\n" +
"  xs = two13 * SQRT(grho)/rho43\n" +
"  xs2 = xs * xs\n" +
"  !\n" +
"  sa2b8 = SQRT(1.0_DP + xs2)\n" +
"  shm1 = LOG(xs + sa2b8)\n" +
"  !\n" +
"  dd = 1.0_DP + 6.0_DP * beta * xs * shm1\n" +
"  dd2 = dd * dd\n" +
"  !\n" +
"  ee = 6.0_DP * beta * xs2 / sa2b8 - 1._DP\n" +
"  sx = two13 * grho / rho43 * ( - beta / dd)\n" +
"  !\n" +
"  v1x = - (4._DP/3._DP) / two13 * xs2 * beta * rho13 * ee / dd2\n" +
"  v2x = two13 * beta * (ee-dd) / (rho43 * dd2)\n" +
"\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE becke88_\n" +
"\n" +
"\n" +
"\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE ggax( rho, grho, sx, v1x, v2x ) !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Perdew-Wang GGA (PW91), exchange part:\n" +
"  !! J.P. Perdew et al.,PRB 46, 6671 (1992)\n" +
"  !\n" +
"  USE kinds, ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, &\n" +
"              dbs, dls\n" +
"  REAL(DP), PARAMETER :: f1=0.19645_DP, f2=7.7956_DP, f3=0.2743_DP, &\n" +
"                         f4=0.1508_DP,  f5=0.004_DP\n" +
"  REAL(DP), PARAMETER :: fp1=-0.019292021296426_DP, fp2=0.161620459673995_DP\n" +
"                       ! fp1= -3/(16 pi)*(3 pi^2)^(-1/3)\n" +
"                       ! fp2= (1/2)(3 pi^2)**(-1/3)\n" +
"  !\n" +
"  rhom43 = rho**(-4.d0/3.d0)\n" +
"  s  = fp2 * SQRT(grho) * rhom43\n" +
"  s2 = s * s\n" +
"  s3 = s2 * s\n" +
"  s4 = s2 * s2\n" +
"  !\n" +
"  exps  = f4 * EXP( - 100.d0 * s2)\n" +
"  as    = f3 - exps - f5 * s2\n" +
"  sa2b8 = SQRT(1.0d0 + f2 * f2 * s2)\n" +
"  shm1  = LOG(f2 * s + sa2b8)\n" +
"  bs    = 1.d0 + f1 * s * shm1 + f5 * s4\n" +
"  !\n" +
"  das = (200.d0 * exps - 2.d0 * f5) * s\n" +
"  dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.d0 * f5 * s3\n" +
"  dls = (das / as - dbs / bs)\n" +
"  !\n" +
"  sx  = fp1 * grho * rhom43 * as / bs\n" +
"  v1x = - 4.d0 / 3.d0 * sx / rho * (1.d0 + s * dls)\n" +
"  v2x = fp1 * rhom43 * as / bs * (2.d0 + s * dls)\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE ggax\n" +
"!\n" +
"!\n" +
"!---------------------------------------------------------------\n" +
"SUBROUTINE pbex( rho, grho, iflag, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------\n" +
"  !! PBE exchange (without Slater exchange):\n" +
"  !! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)\n" +
"  !! iflag=2  \"revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)\n" +
"  !! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)\n" +
"  !! iflag=4  PBEQ2D: L. Chiodo et al., PRL 108, 126402 (2012)\n" +
"  !! iflag=5  optB88: Klimes et al., J. Phys. Cond. Matter, 22, 022201 (2010)\n" +
"  !! iflag=6  optB86b: Klimes et al., Phys. Rev. B 83, 195131 (2011)\n" +
"  !! iflag=7  ev: Engel and Vosko, PRB 47, 13164 (1991)\n" +
"  !! iflag=8  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)\n" +
"  !! iflag=9  W31X: D. Chakraborty, K. Berland, and T. Thonhauser, TBD (2020)\n" +
"  !\n" +
"  USE kinds,      ONLY : DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  INTEGER, INTENT(IN) :: iflag               !<GPU:VALUE>\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  ! input: charge and squared gradient\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  ! output: energy, potential\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx, sx_s\n" +
"  ! (3*pi2*|rho|)^(1/3)\n" +
"  ! |grho|\n" +
"  ! |grho|/(2*kf*|rho|)\n" +
"  ! s^2\n" +
"  ! n*ds/dn\n" +
"  ! n*ds/d(gn)\n" +
"  ! exchange energy LDA part\n" +
"  ! exchange energy gradient part\n" +
"  ! auxiliary variable for energy calculation\n" +
"  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1\n" +
"  REAL(DP) :: p, amu, ab, c, dfxdp, dfxds, upbe, uge, s, ak, aa\n" +
"  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi,        &\n" +
"                         c2=3.093667726280136_DP, c5=4._DP*third, &\n" +
"                         c6=c2*2.51984210_DP, c7=0.8_DP\n" +
"                         ! (3pi^2)^(1/3)*2^(4/3)\n" +
"  ! parameters of the functional\n" +
"  REAL(DP) :: k(9), mu(9), ev(6)\n" +
"  !         pbe         revpbe       pbesol     pbeq2d     optB88   optB86b\n" +
"  !         ev          rpbe         W31x\n" +
"  DATA k  / 0.804_DP,   1.2450_DP,   0.804_DP , 0.804_DP,  1.2_DP,  0.0_DP,       &\n" +
"            0.000_DP,   0.8040_DP,   1.10_DP /,                                   &\n" +
"       mu / 0.2195149727645171_DP, 0.2195149727645171_DP, 0.12345679012345679_DP, &\n" +
"            0.12345679012345679_DP, 0.22_DP, 0.1234_DP, 0.000_DP,                 &\n" +
"            0.2195149727645171_DP, 0.12345679012345679_DP /,                      &\n" +
"       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP,      &\n" +
"            0.011282_DP /  ! a and b parameters of Engel and Vosko\n" +
"  !\n" +
"  SELECT CASE( iflag )\n" +
"  CASE( 4 )\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     p = s1*s1\n" +
"     s = s1\n" +
"     ak = 0.804_DP\n" +
"     amu = 10._DP/81._DP\n" +
"     ab = 0.5217_DP\n" +
"     c = 2._DP\n" +
"     fx =  ak - ak / (1.0_DP + amu * p / ak)  + p**2 * (1.0_DP + p)/       &\n" +
"            (10**c + p**3) * ( -1.0_DP - ak + ak / (1.0_DP + amu * p / ak) &\n" +
"           + ab * p ** (-0.1D1/ 0.4D1) )\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     !\n" +
"     dfxdp = DBLE(1 / (1 + amu * p / ak) ** 2 * amu) + DBLE(2 * p * (1   &\n" +
"     + p) / (10 ** c + p ** 3) * (-1 - ak + ak / (1 + amu * p / ak) + ab &\n" +
"     * p ** (-0.1d1 / 0.4D1))) + DBLE(p ** 2 / (10 ** c + p ** 3) * (    &\n" +
"     -1 - ak + ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) -  &\n" +
"     DBLE(3 * p ** 4 * (1 + p) / (10 ** c + p ** 3) ** 2 * (-1 - ak +    &\n" +
"     ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) + DBLE(p **  &\n" +
"     2) * DBLE(1 + p) / DBLE(10 ** c + p ** 3) * (-DBLE(1 / (1 + amu *   &\n" +
"     p / ak) ** 2 * amu) - DBLE(ab * p ** (-0.5d1 / 0.4D1)) / 0.4D1)\n" +
"     !\n" +
"     dfxds = dfxdp*2._DP*s\n" +
"     dfx = dfxds\n" +
"     ds = - c5 * s1\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  CASE( 5, 9 )\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     ab = mu(iflag) / k(iflag)\n" +
"     p = s1*c6\n" +
"     c = LOG(p + SQRT(p*p+1)) ! asinh(p)\n" +
"     dfx1 = 1 + ab*s1*c\n" +
"     fx = mu(iflag)*s1*s1/dfx1\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     !\n" +
"     dfx = 2*fx/s1-fx/dfx1*(ab*c+ab*s1/SQRT(p*p+1)*c6)\n" +
"     ds  = - c5 * s1\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  CASE( 6 )\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     p = mu(iflag)*s1*s1\n" +
"     fx =  p / ( 1._DP + p )**c7\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     !\n" +
"     dfx = 2*mu(iflag)*s1*fx*(1+(1-c7)*p)/(p*(1+p))\n" +
"     ds = - c5 * s1\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  CASE( 7 )\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     s2 = s1 * s1\n" +
"     s = s2*s2\n" +
"     f1 =  1._DP + ev(1)*s2 + ev(2)*s + ev(3)*s*s2\n" +
"     f2 =  1._DP + ev(4)*s2 + ev(5)*s + ev(6)*s*s2\n" +
"     fx = f1 / f2 - 1._DP\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     ds = - c5 * s1\n" +
"     !\n" +
"     dfx  =  ev(1) + 2*ev(2)*s2 + 3*ev(3)*s\n" +
"     dfx1 =  ev(4) + 2*ev(5)*s2 + 3*ev(6)*s\n" +
"     dfx  = 2 * s1 * ( dfx - f1*dfx1/f2 ) / f2\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  CASE(8)\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     s2 = s1 * s1\n" +
"     f1 = exp( - mu(iflag) * s2 / k(iflag) )\n" +
"     f2 = 1._DP - f1\n" +
"     fx = k(iflag) * f2\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     ds = - c5 * s1\n" +
"     !\n" +
"     dfx = 2._DP * mu(iflag) * s1 * exp( - mu(iflag) * s2 / k(iflag) )\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  CASE DEFAULT\n" +
"     !\n" +
"     agrho = SQRT(grho)\n" +
"     kf = c2 * rho**third\n" +
"     dsg = 0.5_DP / kf\n" +
"     s1 = agrho * dsg / rho\n" +
"     s2 = s1 * s1\n" +
"     f1 = s2 * mu(iflag) / k(iflag)\n" +
"     f2 = 1._DP + f1\n" +
"     f3 = k(iflag) / f2\n" +
"     fx = k(iflag) - f3\n" +
"     !\n" +
"     exunif = - c1 * kf\n" +
"     sx_s = exunif * fx\n" +
"     !\n" +
"     dxunif = exunif * third\n" +
"     ds = - c5 * s1\n" +
"     !\n" +
"     dfx1 = f2 * f2\n" +
"     dfx = 2._DP * mu(iflag) * s1 / dfx1\n" +
"     !\n" +
"     v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"     v2x = exunif * dfx * dsg / agrho\n" +
"     sx  = sx_s * rho\n" +
"     !\n" +
"  END SELECT\n" +
"  !\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE pbex\n" +
"!\n" +
"!\n" +
"!----------------------------------------------------------------------------\n" +
"SUBROUTINE hcth( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !--------------------------------------------------------------------------\n" +
"  !! HCTH/120, JCP 109, p. 6264 (1998)\n" +
"  !! Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)\n" +
"  !! Present release: Mauro Boero, Tsukuba, 11/05/2004\n" +
"  !\n" +
"  !! * rhoa = rhob = 0.5 * rho\n" +
"  !! * grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)\n" +
"  !! * sx  : total exchange correlation energy at point r\n" +
"  !! * v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)\n" +
"  !! * v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)\n" +
"  !\n" +
"  USE kinds,      ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  REAL(DP), PARAMETER :: o3 = 1.0d0/3.0d0, o34 = 4.0d0/3.0d0, fr83 = 8.d0/3.d0\n" +
"  REAL(DP) :: cg0(6), cg1(6), caa(6), cab(6), cx(6)\n" +
"  REAL(DP) :: r3q2, r3pi, gr, rho_o3, rho_o34, xa, xa2, ra, rab,        &\n" +
"       dra_drho, drab_drho, g, dg, era1, dera1_dra, erab0, derab0_drab, &\n" +
"       ex, dex_drho, uaa, uab, ux, ffaa, ffab,  dffaa_drho, dffab_drho, &\n" +
"       denaa, denab, denx, f83rho, bygr, gaa, gab, gx, taa, tab, txx,   &\n" +
"       dgaa_drho, dgab_drho, dgx_drho, dgaa_dgr, dgab_dgr, dgx_dgr\n" +
"  !\n" +
"  r3q2 = 2.d0**(-o3)\n" +
"  r3pi = (3.d0/pi)**o3\n" +
"  ! ... coefficients for pw correlation\n" +
"  cg0(1) = 0.031091d0\n" +
"  cg0(2) = 0.213700d0\n" +
"  cg0(3) = 7.595700d0\n" +
"  cg0(4) = 3.587600d0\n" +
"  cg0(5) = 1.638200d0\n" +
"  cg0(6) = 0.492940d0\n" +
"  cg1(1) = 0.015545d0\n" +
"  cg1(2) = 0.205480d0\n" +
"  cg1(3) =14.118900d0\n" +
"  cg1(4) = 6.197700d0\n" +
"  cg1(5) = 3.366200d0\n" +
"  cg1(6) = 0.625170d0\n" +
"  ! ... hcth-19-4 ...\n" +
"  caa(1) =  0.489508d+00\n" +
"  caa(2) = -0.260699d+00\n" +
"  caa(3) =  0.432917d+00\n" +
"  caa(4) = -0.199247d+01\n" +
"  caa(5) =  0.248531d+01\n" +
"  caa(6) =  0.200000d+00\n" +
"  cab(1) =  0.514730d+00\n" +
"  cab(2) =  0.692982d+01\n" +
"  cab(3) = -0.247073d+02\n" +
"  cab(4) =  0.231098d+02\n" +
"  cab(5) = -0.113234d+02\n" +
"  cab(6) =  0.006000d+00\n" +
"  cx(1)  =  0.109163d+01\n" +
"  cx(2)  = -0.747215d+00\n" +
"  cx(3)  =  0.507833d+01\n" +
"  cx(4)  = -0.410746d+01\n" +
"  cx(5)  =  0.117173d+01\n" +
"  cx(6)  =  0.004000d+00\n" +
"  !  ... ... ... ... ...\n" +
"  !\n" +
"  gr = DSQRT(grho)\n" +
"  rho_o3  = rho**(o3)\n" +
"  rho_o34 = rho**(o34)\n" +
"  xa = 1.25992105d0*gr/rho_o34\n" +
"  xa2 = xa*xa\n" +
"  ra = 0.781592642d0/rho_o3\n" +
"  rab = r3q2*ra\n" +
"  dra_drho = -0.260530881d0/rho_o34\n" +
"  drab_drho = r3q2*dra_drho\n" +
"  CALL pwcorr( ra, cg1, g, dg )                           !<GPU:pwcorr=>pwcorr_d>\n" +
"  era1 = g\n" +
"  dera1_dra = dg\n" +
"  CALL pwcorr( rab, cg0, g, dg )                          !<GPU:pwcorr=>pwcorr_d>\n" +
"  erab0 = g\n" +
"  derab0_drab = dg\n" +
"  ex = -0.75d0*r3pi*rho_o34\n" +
"  dex_drho = -r3pi*rho_o3\n" +
"  uaa = caa(6)*xa2\n" +
"  uaa = uaa/(1.0d0+uaa)\n" +
"  uab = cab(6)*xa2\n" +
"  uab = uab/(1.0d0+uab)\n" +
"  ux = cx(6)*xa2\n" +
"  ux = ux/(1.0d0+ux)\n" +
"  ffaa = rho*era1\n" +
"  ffab = rho*erab0-ffaa\n" +
"  dffaa_drho = era1 + rho*dera1_dra*dra_drho\n" +
"  dffab_drho = erab0 + rho*derab0_drab*drab_drho - dffaa_drho\n" +
"  ! mb-> i-loop removed\n" +
"  denaa = 1.d0 / (1.0d0+caa(6)*xa2)\n" +
"  denab = 1.d0 / (1.0d0+cab(6)*xa2)\n" +
"  denx  = 1.d0 / (1.0d0+cx(6)*xa2)\n" +
"  f83rho = fr83 / rho\n" +
"  bygr = 2.0d0/gr\n" +
"  gaa = caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))\n" +
"  gab = cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))\n" +
"  gx  = cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))\n" +
"  taa = denaa*uaa*(caa(2)+uaa*(2.d0*caa(3)+uaa &\n" +
"        *(3.d0*caa(4)+uaa*4.d0*caa(5))))\n" +
"  tab = denab*uab*(cab(2)+uab*(2.d0*cab(3)+uab &\n" +
"        *(3.d0*cab(4)+uab*4.d0*cab(5))))\n" +
"  txx = denx*ux*(cx(2)+ux*(2.d0*cx(3)+ux &\n" +
"        *(3.d0*cx(4)+ux*4.d0*cx(5))))\n" +
"  dgaa_drho = -f83rho*taa\n" +
"  dgab_drho = -f83rho*tab\n" +
"  dgx_drho  = -f83rho*txx\n" +
"  dgaa_dgr  =  bygr*taa\n" +
"  dgab_dgr  =  bygr*tab\n" +
"  dgx_dgr   =  bygr*txx\n" +
"  ! mb\n" +
"  sx  = ex*gx + ffaa*gaa + ffab*gab\n" +
"  v1x = dex_drho*gx + ex*dgx_drho          &\n" +
"             + dffaa_drho*gaa + ffaa*dgaa_drho &\n" +
"             + dffab_drho*gab + ffab*dgab_drho\n" +
"  v2x = (ex*dgx_dgr + ffaa*dgaa_dgr + ffab*dgab_dgr) / gr\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE hcth\n" +
"    !\n" +
"    !-------------------------------------------------------\n" +
"    SUBROUTINE pwcorr( r, c, g, dg )                    !<GPU:DEVICE>\n" +
"      !-----------------------------------------------------\n" +
"      !\n" +
"      USE kinds,   ONLY: DP\n" +
"      !\n" +
"      IMPLICIT NONE\n" +
"      !\n" +
"      REAL(DP), INTENT(IN)  :: r, c(6)\n" +
"      REAL(DP), INTENT(OUT) :: g, dg\n" +
"      !\n" +
"      ! ... local variables\n" +
"      !\n" +
"      REAL(DP) :: r12, r32, r2, rb, drb, sb\n" +
"      !\n" +
"      r12 = DSQRT(r)\n" +
"      r32 = r*r12\n" +
"      r2  = r*r\n" +
"      rb  = c(3)*r12 + c(4)*r + c(5)*r32 + c(6)*r2\n" +
"      sb  = 1.0d0 + 1.0d0/(2.0d0*c(1)*rb)\n" +
"      g   = -2.0d0 * c(1) * (1.0d0+c(2)*r) * DLOG(sb)\n" +
"      drb = c(3)/(2.0d0*r12) + c(4) + 1.5d0*c(5)*r12 + 2.0d0*c(6)*r\n" +
"      dg  = (1.0d0+c(2)*r)*drb/(rb*rb*sb) - 2.0d0*c(1)*c(2)*DLOG(sb)\n" +
"      !\n" +
"      RETURN\n" +
"      !\n" +
"    END SUBROUTINE pwcorr\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------------\n" +
"SUBROUTINE optx( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------------------\n" +
"  !! OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein\n" +
"  !! Present release: Mauro Boero, Tsukuba, 10/9/2002\n" +
"  !\n" +
"  !! rhoa = rhob = 0.5 * rho in LDA implementation\n" +
"  !! grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)\n" +
"  !! sx  : total exchange correlation energy at point r\n" +
"  !! v1x : d(sx)/drho\n" +
"  !! v2x : 1/gr*d(sx)/d(gr)\n" +
"  !\n" +
"  USE kinds,   ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP), PARAMETER :: small=1.D-30, smal2=1.D-10\n" +
"  ! ... coefficients and exponents\n" +
"  REAL(DP), PARAMETER :: o43=4.0d0/3.0d0, two13=1.259921049894873D0,      &\n" +
"       two53=3.174802103936399D0, gam=0.006D0, a1cx=0.9784571170284421D0, &\n" +
"       a2=1.43169D0\n" +
"  REAL(DP) :: gr, rho43, xa, gamx2, uden, uu\n" +
"  !\n" +
"  ! ... OPTX in compact form\n" +
"  !\n" +
"  gr = MAX(grho,smal2)\n" +
"  rho43 = rho**o43\n" +
"  xa = two13*DSQRT(gr)/rho43\n" +
"  gamx2 = gam*xa*xa\n" +
"  uden = 1.d+00/(1.d+00+gamx2)\n" +
"  uu = a2*gamx2*gamx2*uden*uden\n" +
"  uden = rho43*uu*uden\n" +
"  sx  = -rho43*(a1cx+uu)/two13\n" +
"  v1x = o43*(sx+two53*uden)/rho\n" +
"  v2x = -two53*uden/gr\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE optx\n" +
"!\n" +
"!\n" +
"!---------------------------------------------------------------\n" +
"SUBROUTINE wcx( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------\n" +
"  !!  Wu-Cohen exchange (without Slater exchange):\n" +
"  !!  Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)\n" +
"  !\n" +
"  USE kinds,   ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: kf, agrho, s1, s2, es2, ds, dsg, exunif, fx\n" +
"  ! (3*pi2*|rho|)^(1/3)\n" +
"  ! |grho|\n" +
"  ! |grho|/(2*kf*|rho|)\n" +
"  ! s^2\n" +
"  ! n*ds/dn\n" +
"  ! n*ds/d(gn)\n" +
"  ! exchange energy LDA part\n" +
"  ! exchange energy gradient part\n" +
"  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1, x1, x2, x3, &\n" +
"              dxds1, dxds2, dxds3, sx_s\n" +
"  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  REAL(DP), PARAMETER :: third=1.d0 / 3.d0, c1=0.75d0/pi ,      &\n" +
"                         c2=3.093667726280136d0, c5=4.d0*third, &\n" +
"                         teneightyone = 0.123456790123d0\n" +
"  ! parameters of the functional\n" +
"  REAL(DP), PARAMETER :: k=0.804d0, mu=0.2195149727645171d0, &\n" +
"                         cwc=0.00793746933516d0\n" +
"  !\n" +
"  agrho = SQRT(grho)\n" +
"  kf  = c2 * rho**third\n" +
"  dsg = 0.5d0 / kf\n" +
"  s1  = agrho * dsg / rho\n" +
"  s2  = s1 * s1\n" +
"  es2 = EXP(-s2)\n" +
"  ds  = - c5 * s1\n" +
"  !\n" +
"  !   Energy\n" +
"  ! x = 10/81 s^2 + (mu - 10/81) s^2 e^-s^2 + ln (1 + c s^4)\n" +
"  x1 = teneightyone * s2\n" +
"  x2 = (mu - teneightyone) * s2 * es2\n" +
"  x3 = LOG(1.d0 + cwc * s2 * s2)\n" +
"  f1 = (x1 + x2 + x3) / k\n" +
"  f2 = 1.d0 + f1\n" +
"  f3 = k / f2\n" +
"  fx = k - f3\n" +
"  exunif = - c1 * kf\n" +
"  sx_s = exunif * fx\n" +
"  !\n" +
"  !   Potential\n" +
"  dxunif = exunif * third\n" +
"  dfx1 = f2 * f2\n" +
"  dxds1 = teneightyone\n" +
"  dxds2 = (mu - teneightyone) * es2 * (1.d0 - s2)\n" +
"  dxds3 = 2.d0 * cwc * s2 / (1.d0 + cwc * s2 *s2)\n" +
"  dfx = 2.d0 * s1 * (dxds1 + dxds2 + dxds3) / dfx1\n" +
"  !\n" +
"  v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"  v2x = exunif * dfx * dsg / agrho\n" +
"  sx  = sx_s * rho\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE wcx\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE pbexsr( rho, grho, sxsr, v1xsr, v2xsr, omega )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------------\n" +
"  ! INCLUDE 'cnst.inc'\n" +
"  USE kinds,      ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: omega                !<GPU:VALUE>\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sxsr, v1xsr, v2xsr\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rs, vx, aa, rr, ex, s2, s, d1x, d2x, fx, dsdn, dsdg\n" +
"  REAL(DP), PARAMETER :: small=1.D-20, smal2=1.D-08\n" +
"  REAL(DP), PARAMETER :: us=0.161620459673995492D0, ax=-0.738558766382022406D0, &\n" +
"                         um=0.2195149727645171D0, uk=0.8040D0, ul=um/uk\n" +
"  REAL(DP), PARAMETER :: f1=-1.10783814957303361_DP, alpha=2.0_DP/3.0_DP\n" +
"  !\n" +
"  ! CALL XC(RHO,EX,EC,VX,VC)\n" +
"  !\n" +
"  rs = rho**(1.0_DP/3.0_DP)\n" +
"  vx = (4.0_DP/3.0_DP)*f1*alpha*rs\n" +
"  !\n" +
"  ! aa = dmax1(grho,smal2)\n" +
"  aa = grho\n" +
"  ! rr = rho**(-4.0_DP/3.0_DP)\n" +
"  rr = 1.0_DP/(rho*rs)\n" +
"  ex = ax/rr\n" +
"  s2 = aa*rr*rr*us*us\n" +
"  !\n" +
"  s = SQRT(s2)\n" +
"  IF (s > 8.3D0) THEN\n" +
"     s = 8.572844D0 - 18.796223D0/s2\n" +
"  ENDIF\n" +
"  !\n" +
"  CALL wpbe_analy_erfc_approx_grad( rho, s, omega, fx, d1x, d2x ) !<GPU:wpbe_analy_erfc_approx_grad=>wpbe_analy_erfc_approx_grad_d>\n" +
"  !\n" +
"  sxsr  = ex*fx                        ! - ex\n" +
"  dsdn  = -4.D0/3.D0*s/rho\n" +
"  v1xsr = vx*fx + (dsdn*d2x+d1x)*ex    ! - VX\n" +
"  dsdg  = us*rr\n" +
"  v2xsr = ex*1.D0/SQRT(aa)*dsdg*d2x\n" +
"  !\n" +
"  ! NOTE, here sx is the total energy density,\n" +
"  ! not just the gradient correction energy density as e.g. in pbex()\n" +
"  ! And the same goes for the potentials V1X, V2X\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE pbexsr\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE rPW86( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------------\n" +
"  !! PRB 33, 8800 (1986) and J. Chem. Theory comp. 5, 2754 (2009).\n" +
"  !\n" +
"  USE kinds,      ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, grad_rho, df_ds\n" +
"  REAL(DP), PARAMETER :: a=1.851_DP, b=17.33_DP, c=0.163_DP, &\n" +
"                         s_prefactor=6.18733545256027_DP,    &\n" +
"                         Ax=-0.738558766382022_DP, four_thirds=4._DP/3._DP\n" +
"  !\n" +
"  grad_rho = SQRT(grho)\n" +
"  !\n" +
"  s = grad_rho/(s_prefactor*rho**(four_thirds))\n" +
"  !\n" +
"  s_2 = s**2\n" +
"  s_3 = s_2 * s\n" +
"  s_4 = s_2**2\n" +
"  s_5 = s_3 * s_2\n" +
"  s_6 = s_2 * s_4\n" +
"  !\n" +
"  ! Calculation of energy\n" +
"  fs = (1 + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)\n" +
"  sx = Ax * rho**(four_thirds) * (fs -1.0_DP)\n" +
"  !\n" +
"  ! Calculation of the potential\n" +
"  df_ds = (1._DP/(15._DP*fs**(14.0_DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)\n" +
"  !\n" +
"  v1x = Ax*(four_thirds)*(rho**(1._DP/3._DP)*(fs -1.0_DP) &\n" +
"        -grad_rho/(s_prefactor * rho)*df_ds)\n" +
"  !\n" +
"  v2x = Ax * df_ds/(s_prefactor*grad_rho)\n" +
"  !\n" +
"END SUBROUTINE rPW86\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------\n" +
"SUBROUTINE c09x( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------\n" +
"  !! Cooper '09 exchange for vdW-DF (without Slater exchange):\n" +
"  !! V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)\n" +
"  !\n" +
"  !! Developed thanks to the contribution of\n" +
"  !! Ikutaro Hamada - ikutaro@wpi-aimr.tohoku.ac.jp\n" +
"  !! WPI-Advanced Institute of Materials Research, Tohoku University\n" +
"  !\n" +
"  USE kinds,      ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: kf, agrho, s1, s2, sx_s, ds, dsg, exunif, fx\n" +
"  ! (3*pi2*|rho|)^(1/3)\n" +
"  ! |grho|\n" +
"  ! |grho|/(2*kf*|rho|)\n" +
"  ! s^2\n" +
"  ! n*ds/dn\n" +
"  ! n*ds/d(gn)\n" +
"  ! exchange energy LDA part\n" +
"  ! exchange energy gradient part\n" +
"  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1, dfx2\n" +
"  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi, &\n" +
"                         c2=3.093667726280136_DP, c5=4._DP*third\n" +
"  ! parameters of the functional\n" +
"  REAL(DP) :: kappa, mu, alpha\n" +
"  DATA kappa / 1.245_DP  /, &\n" +
"       mu    / 0.0617_DP /, &\n" +
"       alpha / 0.0483_DP /\n" +
"  !\n" +
"  agrho = SQRT(grho)\n" +
"  kf = c2 * rho**third\n" +
"  dsg = 0.5_DP / kf\n" +
"  s1 = agrho * dsg / rho\n" +
"  s2 = s1 * s1\n" +
"  ds = - c5 * s1\n" +
"  !\n" +
"  ! ... Energy\n" +
"  !\n" +
"  f1 = EXP( - alpha * s2 )\n" +
"  f2 = EXP( - alpha * s2 / 2.0_DP )\n" +
"  f3 = mu * s2 * f1\n" +
"  fx = f3 + kappa * ( 1.0_DP - f2 )\n" +
"  exunif = - c1 * kf\n" +
"  sx_s = exunif * fx\n" +
"  !\n" +
"  ! ... Potential\n" +
"  !\n" +
"  dxunif = exunif * third\n" +
"  dfx1 = 2.0_DP * mu * s1 * ( 1.0_DP - alpha * s2 ) * f1\n" +
"  dfx2 = kappa * alpha * s1 * f2\n" +
"  dfx = dfx1 + dfx2\n" +
"  v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"  v2x = exunif * dfx * dsg / agrho\n" +
"  !\n" +
"  sx  = sx_s * rho\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE c09x\n" +
"!\n" +
"!\n" +
"!---------------------------------------------------------------\n" +
"SUBROUTINE sogga( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-------------------------------------------------------------\n" +
"  !! SOGGA exchange\n" +
"  !\n" +
"  USE kinds,      ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  ! input: charge and abs gradient\n" +
"  ! output: energy and potential\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rho43, xs, xs2, dxs2_drho, dxs2_dgrho2\n" +
"  REAL(DP) :: CX, denom, C1, C2, Fso, Fpbe, ex, Fx, dFx_dxs2, dex_drho\n" +
"  !\n" +
"  REAL(DP), PARAMETER :: one  = 1.0_DP, two   = 2.0_DP, three = 3.0_DP,        &\n" +
"  &                      four = 4.0_DP, eight = 8.0_DP,                        &\n" +
"  &                      f13 = one/three, f23 = two/three,   f43 = four/three, &\n" +
"  &                      f34 = three/four,f83 = eight/three, f12 = one/two\n" +
"  !\n" +
"  REAL(DP), PARAMETER :: mu=0.12346_DP, kapa=0.552_DP\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  !\n" +
"  ! Cx LDA\n" +
"  CX    =  f34 * (three/pi)**f13\n" +
"  denom =  four * (three*pi**two)**f23\n" +
"  C1    =  mu / denom\n" +
"  C2    =  mu / (kapa * denom)\n" +
"  !\n" +
"  rho43 = rho**f43\n" +
"  xs  = grho / rho43\n" +
"  xs2 = xs * xs\n" +
"  !\n" +
"  dxs2_drho   = -f83 * xs2 / rho\n" +
"  dxs2_dgrho2 = one /rho**f83\n" +
"  !\n" +
"  ex       = - CX * rho43\n" +
"  dex_drho = - f43 * CX * rho**f13\n" +
"  !\n" +
"  Fso  = kapa * (one - EXP(-C2*xs2))\n" +
"  Fpbe = C1 * xs2 / (one + C2*xs2)\n" +
"  !\n" +
"  Fx       = f12 * (Fpbe + Fso)\n" +
"  dFx_dxs2 = f12 * (C1 / ((one + C2*xs2)**2) + C1*EXP(-C2*xs2))\n" +
"  !\n" +
"  !   Energy\n" +
"  sx = Fx * ex\n" +
"  !\n" +
"  !   Potential\n" +
"  v1x = dex_drho * Fx  +  ex * dFx_dxs2 * dxs2_drho\n" +
"  v2x = two * ex * dFx_dxs2 * dxs2_dgrho2\n" +
"  !\n" +
"END SUBROUTINE sogga\n" +
"!\n" +
"!\n" +
"!-------------------------------------------------------------------------\n" +
"SUBROUTINE pbexgau( rho, grho, sxsr, v1xsr, v2xsr, alpha_gau )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !\n" +
"  USE kinds,  ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: alpha_gau              !<GPU:VALUE>\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sxsr, v1xsr, v2xsr\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: rs, vx, aa, rr, ex, s2, s, d1x, d2x, fx, dsdn, dsdg\n" +
"  !\n" +
"  REAL(DP), PARAMETER :: small=1.D-20, smal2=1.D-08\n" +
"  REAL(DP), PARAMETER :: us=0.161620459673995492D0, ax=-0.738558766382022406D0, &\n" +
"                         um=0.2195149727645171D0, uk=0.8040D0, ul=um/uk\n" +
"  REAL(DP), PARAMETER :: f1 = -1.10783814957303361_DP, alpha = 2.0_DP/3.0_DP\n" +
"  !\n" +
"  rs = rho**(1.0_DP/3.0_DP)\n" +
"  vx = (4.0_DP/3.0_DP)*f1*alpha*rs\n" +
"  aa = grho\n" +
"  rr = 1.0_DP/(rho*rs)\n" +
"  ex = ax/rr\n" +
"  ! AX is 3/4/PI*(3*PI*PI)**(1/3). This is the same as -c1*c2 in pbex().\n" +
"  s2 = aa*rr*rr*us*us\n" +
"  s = SQRT(s2)\n" +
"  IF (s > 10.D0) THEN\n" +
"     s = 10.D0\n" +
"  ENDIF\n" +
"  CALL pbe_gauscheme( rho, s, alpha_gau, fx, d1x, d2x )   !<GPU:pbe_gauscheme=>pbe_gauscheme_d>\n" +
"  sxsr = ex*fx                        ! - EX\n" +
"  dsdn = -4.D0/3.D0*s/rho\n" +
"  v1xsr = vx*fx + (dsdn*d2x+d1x)*ex   ! - VX\n" +
"  dsdg = us*rr\n" +
"  v2xsr = ex*1.D0/SQRT(aa)*dsdg*d2x\n" +
"  !\n" +
"  ! NOTE, here sx is the total energy density,\n" +
"  ! not just the gradient correction energy density as e.g. in pbex()\n" +
"  ! And the same goes for the potentials V1X, V2X\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE pbexgau\n" +
"    !\n" +
"    !-----------------------------------------------------------------------\n" +
"SUBROUTINE pbe_gauscheme( rho, s, alpha_gau, Fx, dFxdr, dFxds )                    !<GPU:DEVICE>\n" +
"       !--------------------------------------------------------------------\n" +
"       !\n" +
"       IMPLICIT NONE\n" +
"       !\n" +
"       REAL*8 rho,s,alpha_gau,Fx,dFxdr,dFxds\n" +
"       ! input: charge and squared gradient and alpha_gau\n" +
"       ! output: GGA enhancement factor of gau-PBE\n" +
"       ! output: d(Fx)/d(s), d(Fx)/d(rho)\n" +
"       !\n" +
"       REAL*8 Kx, Nx\n" +
"       ! PBE96 GGA enhancement factor\n" +
"       ! GGA enhancement factor of Gaussian Function\n" +
"       !\n" +
"       REAL*8 bx, cx, PI, sqrtpial, Prefac, term_PBE, Third, KsF\n" +
"       REAL*8 d1sdr, d1Kxds, d1Kxdr, d1bxdr, d1bxds, d1bxdKx, &\n" +
"              d1Nxdbx,d1Nxdr, d1Nxds\n" +
"       !\n" +
"       REAL*8 Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten\n" +
"       !\n" +
"       SAVE Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten\n" +
"       DATA Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten &\n" +
"         / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /\n" +
"       !\n" +
"       REAL*8 k , mu\n" +
"       DATA k / 0.804d0 / , mu / 0.21951d0 /\n" +
"       ! parameters of PBE functional\n" +
"       !\n" +
"       Third = One/Three\n" +
"       PI = ACOS(-One)\n" +
"       KsF = (Three*PI*PI*rho)**Third\n" +
"       sqrtpial = SQRT(PI/alpha_gau)\n" +
"       Prefac = Two * SQRT(PI/alpha_gau) / Three\n" +
"       !\n" +
"       ! PBE96 GGA enhancement factor part\n" +
"       term_PBE = One / (One + s*s*mu/k)\n" +
"       Kx =  One + k - k * term_PBE\n" +
"       !\n" +
"       ! GGA enhancement factor of Gaussian Function part\n" +
"       bx = SQRT(Kx*alpha_gau) / KsF\n" +
"       !\n" +
"       ! cx = exp(-One/Four/bx/bx) - One\n" +
"       IF (ABS(One/bx/bx) < 1.0D-4) THEN\n" +
"          cx = TayExp(-One/bx/bx)                               !<GPU:TayExp=>TayExp_d>\n" +
"       ELSE\n" +
"          cx = EXP(-One/bx/bx) - One\n" +
"       ENDIF\n" +
"       !\n" +
"       Nx = bx * Prefac * ( SQRT(PI) * qe_erf(One/bx) + &       !<GPU:qe_erf=>qe_erf_d>\n" +
"        (bx - Two*bx*bx*bx)*cx - Two*bx )\n" +
"       !\n" +
"       ! for convergency\n" +
"       IF (ABS(Nx) < 1.0D-15) THEN\n" +
"         Nx = Zero\n" +
"       ELSEIF ((One - ABS(Nx)) < 1.0D-15) THEN\n" +
"         Nx = One\n" +
"       ELSE\n" +
"         Nx = Nx\n" +
"       ENDIF\n" +
"       ! for convergency end\n" +
"       !\n" +
"       Fx =  Kx * Nx\n" +
"       !\n" +
"       ! 1st derivatives\n" +
"       d1sdr = - Four / Three * s / rho\n" +
"       !\n" +
"       d1Kxds = Two * s * mu * term_PBE * term_PBE\n" +
"       d1Kxdr = d1Kxds * d1sdr\n" +
"       d1bxdKx = bx / (Two* Kx)\n" +
"       !\n" +
"       d1bxdr = - bx /(Three*rho) + d1Kxdr * d1bxdKx\n" +
"       !\n" +
"       d1bxds =  d1bxdKx * d1Kxds\n" +
"       !\n" +
"       d1Nxdbx =  Nx/bx - Prefac * bx * Three * &\n" +
"                   ( cx*(One + Two*bx*bx) + Two )\n" +
"       !\n" +
"       d1Nxdr = d1Nxdbx * d1bxdr\n" +
"       d1Nxds = d1Nxdbx * d1bxds\n" +
"       !\n" +
"       dFxdr = d1Kxdr * Nx + Kx * d1Nxdr\n" +
"       dFxds = d1Kxds * Nx + Kx * d1Nxds\n" +
"       !\n" +
"       RETURN\n" +
"       !\n" +
"END SUBROUTINE pbe_gauscheme\n" +
"!\n" +
"!\n" +
"!-------------------------------------------------\n" +
"FUNCTION TayExp(X)                         !<GPU:DEVICE>\n" +
"  !-------------------------------------------\n" +
"  USE kinds,   ONLY: DP\n" +
"  IMPLICIT NONE\n" +
"  REAL(DP), INTENT(IN) :: X\n" +
"  REAL(DP) :: TAYEXP                        !<GPU:TAYEXP=>TAYEXP_d>\n" +
"  INTEGER :: NTERM,I\n" +
"  REAL(DP) :: SUMVAL,IVAL,COEF\n" +
"  PARAMETER (NTERM=16)\n" +
"  !\n" +
"  SUMVAL = X\n" +
"  IVAL = X\n" +
"  COEF = 1.0D0\n" +
"  DO 10 I = 2, NTERM\n" +
"     COEF = COEF * I\n" +
"     IVAL = IVAL * (X / COEF)\n" +
"     SUMVAL = SUMVAL + IVAL\n" +
"10     CONTINUE\n" +
"  TAYEXP = SUMVAL                      !<GPU:TAYEXP=>TAYEXP_d>\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END FUNCTION TayExp\n" +
"!\n" +
"!\n" +
"!\n" +
"!-------------------------------------------------------------------------\n" +
"SUBROUTINE PW86( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Perdew-Wang 1986 exchange gradient correction: PRB 33, 8800 (1986)\n" +
"  !\n" +
"  USE kinds,  ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, grad_rho, df_ds\n" +
"  REAL(DP), PARAMETER :: a=1.296_DP, b=14._DP, c=0.2_DP,   &\n" +
"                         s_prefactor=6.18733545256027_DP, &\n" +
"                         Ax=-0.738558766382022_DP, four_thirds=4._DP/3._DP\n" +
"  !\n" +
"  grad_rho = SQRT(grho)\n" +
"  !\n" +
"  s = grad_rho / ( s_prefactor*rho**(four_thirds) )\n" +
"  !\n" +
"  s_2 = s**2\n" +
"  s_3 = s_2 * s\n" +
"  s_4 = s_2**2\n" +
"  s_5 = s_3 * s_2\n" +
"  s_6 = s_2 * s_4\n" +
"  !\n" +
"  ! Calculation of energy\n" +
"  fs = (1 + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)\n" +
"  sx = Ax * rho**(four_thirds) * (fs-1._DP)\n" +
"  !\n" +
"  ! Calculation of the potential\n" +
"  df_ds = (1._DP/(15._DP*fs**(14._DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)\n" +
"  !\n" +
"  v1x = Ax*(four_thirds)*( rho**(1._DP/3._DP)*(fs-1._DP) &\n" +
"            -grad_rho/(s_prefactor * rho)*df_ds )\n" +
"  !\n" +
"  v2x = Ax * df_ds/(s_prefactor*grad_rho)\n" +
"  !\n" +
"END SUBROUTINE PW86\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE becke86b( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Becke 1986 gradient correction to exchange\n" +
"  !! A.D. Becke, J. Chem. Phys. 85 (1986) 7184\n" +
"  !\n" +
"  USE kinds, ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: arho, agrho\n" +
"  REAL(DP) :: sgp1, sgp1_45, sgp1_95\n" +
"  REAL(DP) :: rdg2_43, rdg2_73, rdg2_83, rdg2_4, rdg4_5\n" +
"  REAL(DP), PARAMETER :: beta=0.00375_DP, gamma=0.007_DP\n" +
"  !\n" +
"  arho  = 0.5_DP  * rho\n" +
"  agrho = 0.25_DP * grho\n" +
"  !\n" +
"  rdg2_43 = agrho / arho**(4d0/3d0)\n" +
"  rdg2_73 = rdg2_43 / arho\n" +
"  rdg2_83 = rdg2_43 * rdg2_43 / agrho\n" +
"  rdg2_4 = rdg2_43 * rdg2_83 / agrho\n" +
"  rdg4_5 = rdg2_73 * rdg2_83\n" +
"  !\n" +
"  sgp1 = 1d0 + gamma * rdg2_83\n" +
"  sgp1_45 = sgp1**(-4d0/5d0)\n" +
"  sgp1_95 = sgp1_45 / sgp1\n" +
"  !\n" +
"  sx  = -2d0 * beta * agrho / arho**(4d0/3d0) * sgp1_45\n" +
"  v1x = -beta * (-4d0/3d0*rdg2_73*sgp1_45 + 32d0/15d0*gamma*rdg4_5*sgp1_95)\n" +
"  v2x = -beta * (sgp1_45*rdg2_43/agrho - 4d0/5d0 *gamma*rdg2_4*sgp1_95)\n" +
"  !\n" +
"END SUBROUTINE becke86b\n" +
"!\n" +
"!\n" +
"!---------------------------------------------------------------\n" +
"SUBROUTINE b86b( rho, grho, iflag, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-------------------------------------------------------------\n" +
"  !! Becke exchange (without Slater exchange):\n" +
"  !! iflag=1: A. D. Becke, J. Chem. Phys. 85, 7184 (1986) (B86b)\n" +
"  !! iflag=2: J. Klimes, Phys. Rev. B 83, 195131 (2011). (OptB86b)\n" +
"  !! iflag=3: I. Hamada, Phys. Rev. B 89, 121103(R) (B86R)\n" +
"  !! iflag=4: D. Chakraborty, K. Berland, and T. Thonhauser, TBD (2020)\n" +
"  !\n" +
"  !! Ikutaro Hamada - HAMADA.Ikutaro@nims.go.jp\n" +
"  !! National Institute for Materials Science\n" +
"  !\n" +
"  USE kinds,     ONLY : DP\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  INTEGER, INTENT(IN) :: iflag                  !<GPU:VALUE>\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: kf, agrho, s1, s2, sx_s, ds, dsg, exunif, fx\n" +
"  ! (3*pi2*|rho|)^(1/3)\n" +
"  ! |grho|\n" +
"  ! |grho|/(2*kf*|rho|)\n" +
"  ! s^2\n" +
"  ! n*ds/dn\n" +
"  ! n*ds/d(gn)\n" +
"  ! exchange energy LDA part\n" +
"  ! exchange energy gradient part\n" +
"  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1\n" +
"  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )\n" +
"  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0\n" +
"  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi, &\n" +
"                         c2=3.093667726280136_DP, c5=4._DP*third\n" +
"  ! parameters of the functional\n" +
"  REAL(DP) :: k(4), mu(4)\n" +
"  DATA k / 0.5757_DP, 1.0000_DP, 0.711357_DP, 0.58_DP /, &\n" +
"       mu/ 0.2449_DP, 0.1234_DP, 0.1234_DP, 0.12345679012345679_DP /\n" +
"  !\n" +
"  agrho = SQRT(grho)\n" +
"  kf = c2 * rho**third\n" +
"  dsg = 0.5_DP / kf\n" +
"  s1 = agrho * dsg / rho\n" +
"  s2 = s1 * s1\n" +
"  ds = - c5 * s1\n" +
"  !\n" +
"  ! ... Energy\n" +
"  !\n" +
"  f1 = mu(iflag)*s2\n" +
"  f2 = 1._DP + mu(iflag)*s2/k(iflag)\n" +
"  f3 = f2**(4._DP/5._DP)\n" +
"  fx = f1/f3\n" +
"  exunif = - c1 * kf\n" +
"  sx_s = exunif * fx\n" +
"  !\n" +
"  ! ... Potential\n" +
"  !\n" +
"  dxunif = exunif * third\n" +
"  dfx1 = 1._DP + (1._DP/5._DP)*mu(iflag)*s2 / k(iflag)\n" +
"  dfx  = 2._DP * mu(iflag) * s1 * dfx1 / (f2 * f3)\n" +
"  v1x = sx_s + dxunif * fx + exunif * dfx * ds\n" +
"  v2x = exunif * dfx * dsg / agrho\n" +
"  sx = sx_s * rho\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE b86b\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE cx13( rho, grho, sx, v1x, v2x )                    !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! The new exchange partner for a vdW-DF1-cx suggested\n" +
"  !! by K. Berland and P. Hyldgaard, see PRB 89, 035412 (2014),\n" +
"  !! to test the plasmon nature of the vdW-DF1 inner functional.\n" +
"  !\n" +
"  USE kinds, ONLY : DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho, grho\n" +
"  REAL(DP), INTENT(OUT) :: sx, v1x, v2x\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, fs_rPW86, df_rPW86_ds, grad_rho, df_ds\n" +
"  REAL(DP), PARAMETER :: alp=0.021789_DP, beta=1.15_DP, a=1.851_DP, b=17.33_DP, &\n" +
"                         c=0.163_DP, mu_LM=0.09434_DP,    &\n" +
"                         s_prefactor=6.18733545256027_DP, &\n" +
"                         Ax = -0.738558766382022_DP, four_thirds = 4._DP/3._DP\n" +
"  !\n" +
"  grad_rho = SQRT(grho)\n" +
"  !\n" +
"  s = grad_rho/(s_prefactor*rho**(four_thirds))\n" +
"  !\n" +
"  s_2 = s*s\n" +
"  s_3 = s_2 * s\n" +
"  s_4 = s_2 * s_2\n" +
"  s_5 = s_3 * s_2\n" +
"  s_6 = s_2 * s_2 *s_2\n" +
"  !\n" +
"  ! ... Energy\n" +
"  fs_rPW86 = (1._DP + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)\n" +
"  fs = 1._DP/(1._DP + alp*s_6) * (1._DP + mu_LM *s_2) &\n" +
"       + alp*s_6/(beta+alp*s_6)*fs_rPW86\n" +
"  !\n" +
"  sx = Ax * rho**(four_thirds) * (fs-1._DP)\n" +
"  !\n" +
"  ! ... Potential\n" +
"  df_rPW86_ds = (1._DP/(15._DP*fs_rPW86**(14._DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)\n" +
"  !\n" +
"  df_ds = 1._DP/(1._DP+alp*s_6)**2*( 2._DP*mu_LM*s*(1._DP+alp*s_6) &\n" +
"            - 6._DP*alp*s_5*( 1._DP+mu_LM*s_2) )                   &\n" +
"          + alp*s_6/(beta+alp*s_6)*df_rPW86_ds                     &\n" +
"          + 6._DP*alp*s_5*fs_rPW86/(beta+alp*s_6)*(1._DP-alp*s_6/(beta + alp*s_6))\n" +
"  !\n" +
"  v1x = Ax*(four_thirds)*(rho**(1._DP/3._DP)*(fs-1._DP) &\n" +
"        -grad_rho/(s_prefactor * rho)*df_ds)\n" +
"  v2x = Ax * df_ds/(s_prefactor*grad_rho)\n" +
"  !\n" +
"END SUBROUTINE cx13\n" +
"!\n" +
"!\n" +
"!\n" +
"! ===========> SPIN <===========\n" +
"!\n" +
"!-----------------------------------------------------------------------\n" +
"SUBROUTINE becke88_spin_( rho_up, rho_dw, grho_up, grho_dw, sx_up, sx_dw, v1x_up, v1x_dw, v2x_up, v2x_dw )                     !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case\n" +
"  !\n" +
"  USE kinds,    ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho_up, rho_dw\n" +
"  !! charge\n" +
"  REAL(DP), INTENT(IN) :: grho_up, grho_dw\n" +
"  !! gradient\n" +
"  REAL(DP), INTENT(OUT) :: sx_up, sx_dw\n" +
"  !! the up and down energies\n" +
"  REAL(DP), INTENT(OUT) :: v1x_up, v1x_dw\n" +
"  !! first part of the potential\n" +
"  REAL(DP), INTENT(OUT) :: v2x_up, v2x_dw\n" +
"  !! second part of the potential\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  INTEGER :: is\n" +
"  REAL(DP), PARAMETER :: beta = 0.0042_DP, third = 1._DP/3._DP\n" +
"  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee\n" +
"  !\n" +
"  !\n" +
"  !DO is = 1, 2\n" +
"     CALL ann_b88(rho_up, grho_up, sx_up, v1x_up, v2x_up)\n" +
"     CALL ann_b88(rho_dw, grho_dw, sx_dw, v1x_dw, v2x_dw)\n" +
"\n" +
"  !ENDDO\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE becke88_spin_\n" +
"\n" +
"SUBROUTINE becke88_spin( rho_up, rho_dw, grho_up, grho_dw, sx_up, sx_dw, v1x_up, v1x_dw, v2x_up, v2x_dw )                     !<GPU:DEVICE>\n" +
"  !-----------------------------------------------------------------------\n" +
"  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case\n" +
"  !\n" +
"  USE kinds,    ONLY: DP\n" +
"  !\n" +
"  IMPLICIT NONE\n" +
"  !\n" +
"  REAL(DP), INTENT(IN) :: rho_up, rho_dw\n" +
"  !! charge\n" +
"  REAL(DP), INTENT(IN) :: grho_up, grho_dw\n" +
"  !! gradient\n" +
"  REAL(DP), INTENT(OUT) :: sx_up, sx_dw\n" +
"  !! the up and down energies\n" +
"  REAL(DP), INTENT(OUT) :: v1x_up, v1x_dw\n" +
"  !! first part of the potential\n" +
"  REAL(DP), INTENT(OUT) :: v2x_up, v2x_dw\n" +
"  !! second part of the potential\n" +
"  !\n" +
"  ! ... local variables\n" +
"  !\n" +
"  INTEGER :: is\n" +
"  REAL(DP), PARAMETER :: beta = 0.0042_DP, third = 1._DP/3._DP\n" +
"  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee\n" +
"  !\n" +
"  !\n" +
"\n" +
"  !DO is = 1, 2\n" +
"     rho13 = rho_up**third\n" +
"     rho43 = rho13**4\n" +
"     xs  = SQRT(grho_up) / rho43\n" +
"     xs2 = xs * xs\n" +
"     sa2b8 = SQRT(1.0d0 + xs2)\n" +
"     shm1  = LOG(xs + sa2b8)\n" +
"     dd  = 1.0d0 + 6.0d0 * beta * xs * shm1\n" +
"     dd2 = dd * dd\n" +
"     ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0\n" +
"     sx_up  = grho_up / rho43 * (-beta/dd)\n" +
"     v1x_up = -(4.d0/3.d0) * xs2 * beta * rho13 * ee / dd2\n" +
"     v2x_up = beta * (ee-dd) / (rho43*dd2)\n" +
"\n" +
"     rho13 = rho_dw**third\n" +
"     rho43 = rho13**4\n" +
"     xs  = SQRT(grho_dw) / rho43\n" +
"     xs2 = xs * xs\n" +
"     sa2b8 = SQRT(1.0d0 + xs2)\n" +
"     shm1  = LOG(xs + sa2b8)\n" +
"     dd  = 1.0d0 + 6.0d0 * beta * xs * shm1\n" +
"     dd2 = dd * dd\n" +
"     ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0\n" +
"     sx_dw  = grho_dw / rho43 * (-beta/dd)\n" +
"     v1x_dw = -(4.d0/3.d0) * xs2 * beta * rho13 * ee / dd2\n" +
"     v2x_dw = beta * (ee-dd) / (rho43*dd2)\n" +
"  !ENDDO\n" +
"  !\n" +
"  RETURN\n" +
"  !\n" +
"END SUBROUTINE becke88_spin\n" +
"!\n" +
"!\n" +
"!-----------------------------------------------------------------------------\n" +
"SUBROUTINE wpbe_analy_erfc_approx_grad( rho, s, omega, Fx_wpbe, d1rfx, d1sfx )                     !<GPU:DEVICE>\n" +
"      !-----------------------------------------------------------------------\n" +
"      !\n" +
"      !     wPBE Enhancement Factor (erfc approx.,analytical, gradients)\n" +
"      !\n" +
"      !--------------------------------------------------------------------\n" +
"      !\n" +
"      USE kinds,    ONLY: DP\n" +
"      IMPLICIT NONE\n" +
"      !\n" +
"      REAL(DP) rho,s,omega,Fx_wpbe,d1sfx,d1rfx\n" +
"      !\n" +
"      REAL(DP) f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516,f98\n" +
"      REAL(DP) pi,pi2,pi_23,srpi\n" +
"      REAL(DP) Three_13\n" +
"      !\n" +
"      REAL(DP) ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8\n" +
"      REAL(DP) eb1\n" +
"      REAL(DP) A,B,C,D,E\n" +
"      REAL(DP) Ha1,Ha2,Ha3,Ha4,Ha5\n" +
"      REAL(DP) Fc1,Fc2\n" +
"      REAL(DP) EGa1,EGa2,EGa3\n" +
"      REAL(DP) EGscut,wcutoff,expfcutoff\n" +
"      !\n" +
"      REAL(DP) xkf, xkfrho\n" +
"      REAL(DP) w,w2,w3,w4,w5,w6,w7,w8\n" +
"      REAL(DP) d1rw\n" +
"      REAL(DP) A2,A3,A4,A12,A32,A52,A72\n" +
"      REAL(DP) X\n" +
"      REAL(DP) s2,s3,s4,s5,s6\n" +
"      !\n" +
"      REAL(DP) H,F\n" +
"      REAL(DP) Hnum,Hden,d1sHnum,d1sHden\n" +
"      REAL(DP) d1sH,d1sF\n" +
"      REAL(DP) G_a,G_b,EG\n" +
"      REAL(DP) d1sG_a,d1sG_b,d1sEG\n" +
"      !\n" +
"      REAL(DP) Hsbw,Hsbw2,Hsbw3,Hsbw4,Hsbw12,Hsbw32,Hsbw52,Hsbw72\n" +
"      REAL(DP) DHsbw,DHsbw2,DHsbw3,DHsbw4,DHsbw5\n" +
"      REAL(DP) DHsbw12,DHsbw32,DHsbw52,DHsbw72,DHsbw92\n" +
"      REAL(DP) d1sHsbw,d1rHsbw\n" +
"      REAL(DP) d1sDHsbw,d1rDHsbw\n" +
"      REAL(DP) HsbwA94,HsbwA9412\n" +
"      REAL(DP) HsbwA942,HsbwA943,HsbwA945\n" +
"      REAL(DP) piexperf,expei\n" +
"      REAL(DP) piexperfd1,expeid1\n" +
"      REAL(DP) d1spiexperf,d1sexpei\n" +
"      REAL(DP) d1rpiexperf,d1rexpei\n" +
"      REAL(DP) expei1,expei2,expei3,expei4\n" +
"      !\n" +
"      REAL(DP) DHs,DHs2,DHs3,DHs4,DHs72,DHs92,DHsw,DHsw2,DHsw52,DHsw72\n" +
"      REAL(DP) d1sDHs,d1rDHsw\n" +
"      !\n" +
"      REAL(DP) np1,np2\n" +
"      REAL(DP) d1rnp1,d1rnp2\n" +
"      REAL(DP) t1,t2t9,t10,t10d1\n" +
"      REAL(DP) f2,f3,f4,f5,f6,f7,f8,f9\n" +
"      REAL(DP) f2d1,f3d1,f4d1,f5d1,f6d1,f8d1,f9d1\n" +
"      REAL(DP) d1sf2,d1sf3,d1sf4,d1sf5,d1sf6,d1sf7,d1sf8,d1sf9\n" +
"      REAL(DP) d1rf2,d1rf3,d1rf4,d1rf5,d1rf6,d1rf7,d1rf8,d1rf9\n" +
"      REAL(DP) d1st1,d1rt1\n" +
"      REAL(DP) d1st2t9,d1rt2t9\n" +
"      REAL(DP) d1st10,d1rt10\n" +
"      REAL(DP) d1sterm1,d1rterm1,term1d1\n" +
"      REAL(DP) d1sterm2\n" +
"      REAL(DP) d1sterm3,d1rterm3\n" +
"      REAL(DP) d1sterm4,d1rterm4\n" +
"      REAL(DP) d1sterm5,d1rterm5\n" +
"      !\n" +
"      REAL(DP) term1,term2,term3,term4,term5\n" +
"      !\n" +
"      REAL(DP) ax,um,uk,ul\n" +
"      REAL(DP) gc1,gc2\n" +
"      !\n" +
"      ! REAL(DP) ei\n" +
"      !\n" +
"      REAL(DP) Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten\n" +
"      REAL(DP) Fifteen,Sixteen\n" +
"      REAL(DP) r12,r64,r36,r81,r256,r384,r864,r1944,r4374\n" +
"      REAL(DP) r20,r25,r27,r48,r120,r128,r144,r288,r324,r512,r729\n" +
"      REAL(DP) r30,r32,r75,r243,r2187,r6561,r40,r105,r54,r135\n" +
"      REAL(DP) r1215,r15309\n" +
"      !\n" +
"      SAVE Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten\n" +
"      DATA Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten &\n" +
"        / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /\n" +
"      SAVE Fifteen,Sixteen\n" +
"      DATA Fifteen,Sixteen / 1.5D1, 1.6D1 /\n" +
"      SAVE r36,r64,r81,r256,r384,r864,r1944,r4374\n" +
"      DATA r36,r64,r81,r256,r384,r864,r1944,r4374 &\n" +
"        / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2,8.64D2,1.944D3,4.374D3 /\n" +
"      SAVE r27,r48,r120,r128,r144,r288,r324,r512,r729\n" +
"      DATA r27,r48,r120,r128,r144,r288,r324,r512,r729 &\n" +
"        / 2.7D1,4.8D1,1.2D2,1.28D2,1.44D2,2.88D2,3.24D2,5.12D2,7.29D2 /\n" +
"      SAVE r20,r32,r243,r2187,r6561,r40\n" +
"      DATA r20,r32,r243,r2187,r6561,r40 &\n" +
"        / 2.0d1,3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /\n" +
"      SAVE r12,r25,r30,r54,r75,r105,r135,r1215,r15309\n" +
"      DATA r12,r25,r30,r54,r75,r105,r135,r1215,r15309 &\n" +
"        / 1.2D1,2.5d1,3.0d1,5.4D1,7.5d1,1.05D2,1.35D2,1.215D3,1.5309D4 /\n" +
"      !\n" +
"      ! ... General constants\n" +
"      !\n" +
"      f12    = 0.5d0\n" +
"      f13    = One/Three\n" +
"      f14    = 0.25d0\n" +
"      f18    = 0.125d0\n" +
"      !\n" +
"      f23    = Two * f13\n" +
"      f43    = Two * f23\n" +
"      !\n" +
"      f32    = 1.5d0\n" +
"      f72    = 3.5d0\n" +
"      f34    = 0.75d0\n" +
"      f94    = 2.25d0\n" +
"      f98    = 1.125d0\n" +
"      f1516  = Fifteen / Sixteen\n" +
"      !\n" +
"      pi     = ACOS(-One)\n" +
"      pi2    = pi*pi\n" +
"      pi_23  = pi2**f13\n" +
"      srpi   = SQRT(pi)\n" +
"      !\n" +
"      Three_13 = Three**f13\n" +
"      !\n" +
"      ! Constants from fit\n" +
"      !\n" +
"      ea1 = -1.128223946706117d0\n" +
"      ea2 = 1.452736265762971d0\n" +
"      ea3 = -1.243162299390327d0\n" +
"      ea4 = 0.971824836115601d0\n" +
"      ea5 = -0.568861079687373d0\n" +
"      ea6 = 0.246880514820192d0\n" +
"      ea7 = -0.065032363850763d0\n" +
"      ea8 = 0.008401793031216d0\n" +
"      !\n" +
"      eb1 = 1.455915450052607d0\n" +
"      !\n" +
"      ! Constants for PBE hole\n" +
"      !\n" +
"      A      =  1.0161144d0\n" +
"      B      = -3.7170836d-1\n" +
"      C      = -7.7215461d-2\n" +
"      D      =  5.7786348d-1\n" +
"      E      = -5.1955731d-2\n" +
"      X      = - Eight/Nine\n" +
"      !\n" +
"      ! Constants for fit of H(s) (PBE)\n" +
"      !\n" +
"      Ha1    = 9.79681d-3\n" +
"      Ha2    = 4.10834d-2\n" +
"      Ha3    = 1.87440d-1\n" +
"      Ha4    = 1.20824d-3\n" +
"      Ha5    = 3.47188d-2\n" +
"      !\n" +
"      ! Constants for F(H) (PBE)\n" +
"      !\n" +
"      Fc1    = 6.4753871d0\n" +
"      Fc2    = 4.7965830d-1\n" +
"      !\n" +
"      ! Constants for polynomial expansion for EG for small s\n" +
"      !\n" +
"      EGa1   = -2.628417880d-2\n" +
"      EGa2   = -7.117647788d-2\n" +
"      EGa3   =  8.534541323d-2\n" +
"      !\n" +
"      ! Constants for large x expansion of exp(x)*ei(-x)\n" +
"      !\n" +
"      expei1 = 4.03640D0\n" +
"      expei2 = 1.15198D0\n" +
"      expei3 = 5.03627D0\n" +
"      expei4 = 4.19160D0\n" +
"      !\n" +
"      ! Cutoff criterion below which to use polynomial expansion\n" +
"      !\n" +
"      EGscut     = 8.0d-2\n" +
"      wcutoff    = 1.4D1\n" +
"      expfcutoff = 7.0D2\n" +
"      !\n" +
"      ! Calculate prelim variables\n" +
"      !\n" +
"      xkf    = (Three*pi2*rho) ** f13\n" +
"      xkfrho = xkf * rho\n" +
"      !\n" +
"      A2 = A*A\n" +
"      A3 = A2*A\n" +
"      A4 = A3*A\n" +
"      A12 = SQRT(A)\n" +
"      A32 = A12*A\n" +
"      A52 = A32*A\n" +
"      A72 = A52*A\n" +
"      !\n" +
"      w      = omega / xkf\n" +
"      w2    = w * w\n" +
"      w3    = w2 * w\n" +
"      w4    = w2 * w2\n" +
"      w5    = w3 * w2\n" +
"      w6    = w5 * w\n" +
"      w7    = w6 * w\n" +
"      w8    = w7 * w\n" +
"      !\n" +
"      d1rw  = -(One/(Three*rho))*w\n" +
"      !\n" +
"      X      = - Eight/Nine\n" +
"      !\n" +
"      s2     = s*s\n" +
"      s3     = s2*s\n" +
"      s4     = s2*s2\n" +
"      s5     = s4*s\n" +
"      s6     = s5*s\n" +
"      !\n" +
"      ! Calculate wPBE enhancement factor\n" +
"      !\n" +
"      Hnum    = Ha1*s2 + Ha2*s4\n" +
"      Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6\n" +
"      !\n" +
"      H       = Hnum/Hden\n" +
"      !\n" +
"      d1sHnum = Two*Ha1*s + Four*Ha2*s3\n" +
"      d1sHden = Four*Ha3*s3 + Five*Ha4*s4 + Six*Ha5*s5\n" +
"      !\n" +
"      d1sH    = (Hden*d1sHnum - Hnum*d1sHden) / (Hden*Hden)\n" +
"      !\n" +
"      F      = Fc1*H + Fc2\n" +
"      d1sF   = Fc1*d1sH\n" +
"      !\n" +
"      ! Change exponent of Gaussian if we're using the simple approx.\n" +
"      !\n" +
"      IF (w > wcutoff) eb1 = 2.0d0\n" +
"      !\n" +
"      ! Calculate helper variables (should be moved later on...)\n" +
"      !\n" +
"      Hsbw = s2*H + eb1*w2\n" +
"      Hsbw2 = Hsbw*Hsbw\n" +
"      Hsbw3 = Hsbw2*Hsbw\n" +
"      Hsbw4 = Hsbw3*Hsbw\n" +
"      Hsbw12 = SQRT(Hsbw)\n" +
"      Hsbw32 = Hsbw12*Hsbw\n" +
"      Hsbw52 = Hsbw32*Hsbw\n" +
"      Hsbw72 = Hsbw52*Hsbw\n" +
"      !\n" +
"      d1sHsbw  = d1sH*s2 + Two*s*H\n" +
"      d1rHsbw  = Two*eb1*d1rw*w\n" +
"      !\n" +
"      DHsbw = D + s2*H + eb1*w2\n" +
"      DHsbw2 = DHsbw*DHsbw\n" +
"      DHsbw3 = DHsbw2*DHsbw\n" +
"      DHsbw4 = DHsbw3*DHsbw\n" +
"      DHsbw5 = DHsbw4*DHsbw\n" +
"      DHsbw12 = SQRT(DHsbw)\n" +
"      DHsbw32 = DHsbw12*DHsbw\n" +
"      DHsbw52 = DHsbw32*DHsbw\n" +
"      DHsbw72 = DHsbw52*DHsbw\n" +
"      DHsbw92 = DHsbw72*DHsbw\n" +
"      !\n" +
"      HsbwA94   = f94 * Hsbw / A\n" +
"      HsbwA942  = HsbwA94*HsbwA94\n" +
"      HsbwA943  = HsbwA942*HsbwA94\n" +
"      HsbwA945  = HsbwA943*HsbwA942\n" +
"      HsbwA9412 = SQRT(HsbwA94)\n" +
"      !\n" +
"      DHs    = D + s2*H\n" +
"      DHs2   = DHs*DHs\n" +
"      DHs3   = DHs2*DHs\n" +
"      DHs4   = DHs3*DHs\n" +
"      DHs72  = DHs3*SQRT(DHs)\n" +
"      DHs92  = DHs72*DHs\n" +
"      !\n" +
"      d1sDHs = Two*s*H + s2*d1sH\n" +
"      !\n" +
"      DHsw   = DHs + w2\n" +
"      DHsw2  = DHsw*DHsw\n" +
"      DHsw52 = SQRT(DHsw)*DHsw2\n" +
"      DHsw72 = DHsw52*DHsw\n" +
"      !\n" +
"      d1rDHsw = Two*d1rw*w\n" +
"      !\n" +
"      IF (s > EGscut) THEN\n" +
"        !\n" +
"        G_a    = srpi * (Fifteen*E + Six*C*(One+F*s2)*DHs + &\n" +
"                         Four*B*(DHs2) + Eight*A*(DHs3))    &\n" +
"                      * (One / (Sixteen * DHs72))           &\n" +
"                       - f34*pi*SQRT(A) * EXP(f94*H*s2/A) * &\n" +
"                         (One - qe_erf(f32*s*SQRT(H/A)))                    !<GPU:qe_erf=>qe_erf_d>\n" +
"        !\n" +
"        d1sG_a = (One/r32)*srpi *                           &\n" +
"                 ((r36*(Two*H + d1sH*s) / (A12*SQRT(H/A)))  &\n" +
"                  + (One/DHs92) *                           &\n" +
"                     (-Eight*A*d1sDHs*DHs3 - r105*d1sDHs*E  &\n" +
"                      -r30*C*d1sDHs*DHs*(One+s2*F)          &\n" +
"                      +r12*DHs2*(-B*d1sDHs + C*s*(d1sF*s + Two*F)))  &\n" +
"                  - ((r54*EXP(f94*H*s2/A)*srpi*s*(Two*H+d1sH*s)*     &\n" +
"                     qe_erfc(f32*SQRT(H/A)*s))                       &      !<GPU:qe_erfc=>qe_erfc_d>\n" +
"                     / A12))\n" +
"        !\n" +
"        G_b    = (f1516 * srpi * s2) / DHs72\n" +
"        !\n" +
"        d1sG_b = (Fifteen*srpi*s*(Four*DHs - Seven*d1sDHs*s)) &\n" +
"                 / (r32*DHs92)\n" +
"        !\n" +
"        EG     = - (f34*pi + G_a) / G_b\n" +
"        !\n" +
"        d1sEG  = (-Four*d1sG_a*G_b + d1sG_b*(Four*G_a + Three*pi)) &\n" +
"                 / (Four*G_b*G_b)\n" +
"        !\n" +
"      ELSE\n" +
"        !\n" +
"        EG    = EGa1 + EGa2*s2 + EGa3*s4\n" +
"        d1sEG = Two*EGa2*s + Four*EGa3*s3\n" +
"        !\n" +
"      ENDIF\n" +
"      !\n" +
"      ! Calculate the terms needed in any case\n" +
"      !\n" +
"      term2 =       (DHs2*B + DHs*C + Two*E + DHs*s2*C*F + Two*s2*EG) / &\n" +
"                    (Two*DHs3)\n" +
"      !\n" +
"      d1sterm2 = (-Six*d1sDHs*(EG*s2 + E)                     &\n" +
"                  + DHs2 * (-d1sDHs*B + s*C*(d1sF*s + Two*F)) &\n" +
"                  + Two*DHs * (Two*EG*s - d1sDHs*C            &\n" +
"                  + s2 * (d1sEG - d1sDHs*C*F)))               &\n" +
"                 / (Two*DHs4)\n" +
"\n" +
"      term3 = - w  * (Four*DHsw2*B + Six*DHsw*C + Fifteen*E &\n" +
"                      + Six*DHsw*s2*C*F + Fifteen*s2*EG) /  &\n" +
"                     (Eight*DHs*DHsw52)\n" +
"      !\n" +
"      d1sterm3 = w * (Two*d1sDHs*DHsw * (Four*DHsw2*B         &\n" +
"                         + Six*DHsw*C + Fifteen*E             &\n" +
"                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &\n" +
"                      + DHs * (r75*d1sDHs*(EG*s2 + E)         &\n" +
"                         + Four*DHsw2*(d1sDHs*B               &\n" +
"                              - Three*s*C*(d1sF*s + Two*F))   &\n" +
"                         - Six*DHsw*(-Three*d1sDHs*C          &\n" +
"                              + s*(Ten*EG + Five*d1sEG*s      &\n" +
"                                  - Three*d1sDHs*s*C*F))))    &\n" +
"                 / (Sixteen*DHs2*DHsw72)\n" +
"      !\n" +
"      d1rterm3 = (-Two*d1rw*DHsw * (Four*DHsw2*B              &\n" +
"                         + Six*DHsw*C + Fifteen*E             &\n" +
"                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &\n" +
"                      + w * d1rDHsw * (r75*(EG*s2 + E)        &\n" +
"                         + Two*DHsw*(Two*DHsw*B + Nine*C      &\n" +
"                                     + Nine*s2*C*F)))         &\n" +
"                 / (Sixteen*DHs*DHsw72)\n" +
"\n" +
"      term4 = - w3 * (DHsw*C + Five*E + DHsw*s2*C*F + Five*s2*EG) /  &\n" +
"                     (Two*DHs2*DHsw52)\n" +
"      !\n" +
"      d1sterm4 = (w3 * (Four*d1sDHs*DHsw * (DHsw*C + Five*E   &\n" +
"                             + s2 * (Five*EG + DHsw*C*F))     &\n" +
"                        + DHs * (r25*d1sDHs*(EG*s2 + E)       &\n" +
"                             - Two*DHsw2*s*C*(d1sF*s + Two*F) &\n" +
"                             + DHsw * (Three*d1sDHs*C + s*(-r20*EG  &\n" +
"                                   - Ten*d1sEG*s              &\n" +
"                                   + Three*d1sDHs*s*C*F)))))  &\n" +
"                 / (Four*DHs3*DHsw72)\n" +
"      !\n" +
"      d1rterm4 = (w2 * (-Six*d1rw*DHsw * (DHsw*C + Five*E   &\n" +
"                             + s2 * (Five*EG + DHsw*C*F))   &\n" +
"                        + w * d1rDHsw * (r25*(EG*s2 + E) +  &\n" +
"                             Three*DHsw*C*(One + s2*F))))  &\n" +
"                 / (Four*DHs2*DHsw72)\n" +
"      !\n" +
"      term5 = - w5 * (E + s2*EG) / &\n" +
"                     (DHs3*DHsw52)\n" +
"      !\n" +
"      d1sterm5 = (w5 * (Six*d1sDHs*DHsw*(EG*s2 + E)               &\n" +
"                        + DHs * (-Two*DHsw*s * (Two*EG + d1sEG*s) &\n" +
"                             + Five*d1sDHs * (EG*s2 + E))))       &\n" +
"                 / (Two*DHs4*DHsw72)\n" +
"      !\n" +
"      d1rterm5 = (w4 * Five*(EG*s2 + E) * (-Two*d1rw*DHsw   &\n" +
"                                           + d1rDHsw * w))  &\n" +
"                 / (Two*DHs3*DHsw72)\n" +
"      !\n" +
"      !\n" +
"      IF ((s > 0.0d0).OR.(w > 0.0d0)) THEN\n" +
"        !\n" +
"        t10    = (f12)*A*LOG(Hsbw / DHsbw)\n" +
"        t10d1  = f12*A*(One/Hsbw - One/DHsbw)\n" +
"        d1st10 = d1sHsbw*t10d1\n" +
"        d1rt10 = d1rHsbw*t10d1\n" +
"        !\n" +
"      ENDIF\n" +
"      !\n" +
"      ! Calculate exp(x)*f(x) depending on size of x\n" +
"      !\n" +
"      IF (HsbwA94 < expfcutoff) THEN\n" +
"        !\n" +
"        piexperf = pi*EXP(HsbwA94)*qe_erfc(HsbwA9412)                   !<GPU:qe_erfc=>qe_erfc_d>\n" +
"        ! expei    = Exp(HsbwA94)*Ei(-HsbwA94)\n" +
"        expei    = EXP(HsbwA94)*(-expint(1,HsbwA94))                   !<GPU:expint=>expint_d>\n" +
"\n" +
"      ELSE\n" +
"        !\n" +
"        ! print *,rho,s,\" LARGE HsbwA94\"\n" +
"        !\n" +
"        piexperf = pi*(One/(srpi*HsbwA9412)          &\n" +
"                   - One/(Two*SQRT(pi*HsbwA943))     &\n" +
"                   + Three/(Four*SQRT(pi*HsbwA945)))\n" +
"        !\n" +
"        expei  = - (One/HsbwA94) *                         &\n" +
"                   (HsbwA942 + expei1*HsbwA94 + expei2) /  &\n" +
"                   (HsbwA942 + expei3*HsbwA94 + expei4)\n" +
"\n" +
"      ENDIF\n" +
"      !\n" +
"      ! Calculate the derivatives (based on the orig. expression)\n" +
"      ! --> Is this ok? ==> seems to be ok...\n" +
"      !\n" +
"      piexperfd1  = - (Three*srpi*SQRT(Hsbw/A))/(Two*Hsbw)  &\n" +
"                    + (Nine*piexperf)/(Four*A)\n" +
"      d1spiexperf = d1sHsbw*piexperfd1\n" +
"      d1rpiexperf = d1rHsbw*piexperfd1\n" +
"\n" +
"      expeid1  = f14*(Four/Hsbw + (Nine*expei)/A)\n" +
"      d1sexpei = d1sHsbw*expeid1\n" +
"      d1rexpei = d1rHsbw*expeid1\n" +
"      !\n" +
"      IF (w == Zero) THEN\n" +
"        !\n" +
"        ! Fall back to original expression for the PBE hole\n" +
"        !\n" +
"        t1 = -f12*A*expei\n" +
"        d1st1 = -f12*A*d1sexpei\n" +
"        d1rt1 = -f12*A*d1rexpei\n" +
"        !\n" +
"        ! write(*,*) s, t1, t10, d1st1,d1rt1,d1rt10\n" +
"        !\n" +
"        IF (s > 0.0D0) THEN\n" +
"          !\n" +
"          term1    = t1 + t10\n" +
"          d1sterm1 = d1st1 + d1st10\n" +
"          d1rterm1 = d1rt1 + d1rt10\n" +
"          !\n" +
"          Fx_wpbe = X * (term1 + term2)\n" +
"          !\n" +
"          d1sfx = X * (d1sterm1 + d1sterm2)\n" +
"          d1rfx = X * d1rterm1\n" +
"          !\n" +
"        ELSE\n" +
"          !\n" +
"          Fx_wpbe = 1.0d0\n" +
"          !\n" +
"          ! TODO    This is checked to be true for term1\n" +
"          !         How about the other terms???\n" +
"          !\n" +
"          d1sfx   = 0.0d0\n" +
"          d1rfx   = 0.0d0\n" +
"          !\n" +
"        ENDIF\n" +
"        !\n" +
"        !\n" +
"      ELSEIF (w > wcutoff) THEN\n" +
"        !\n" +
"        ! Use simple Gaussian approximation for large w\n" +
"        !\n" +
"        ! print *,rho,s,\" LARGE w\"\n" +
"        !\n" +
"        term1 = -f12*A*(expei+LOG(DHsbw)-LOG(Hsbw))\n" +
"\n" +
"        term1d1  = - A/(Two*DHsbw) - f98*expei\n" +
"        d1sterm1 = d1sHsbw*term1d1\n" +
"        d1rterm1 = d1rHsbw*term1d1\n" +
"\n" +
"        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)\n" +
"\n" +
"        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3  &\n" +
"                              + d1sterm4 + d1sterm5)\n" +
"\n" +
"        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)\n" +
"        !\n" +
"      ELSE\n" +
"         !\n" +
"         ! For everything else, use the full blown expression\n" +
"         !\n" +
"         ! First, we calculate the polynomials for the first term\n" +
"         !\n" +
"         np1    = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)     &\n" +
"                  - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)\n" +
"        !\n" +
"        d1rnp1 = - f32*ea1*d1rw*A12 + (r81*ea3*d1rw*w2)/(Eight*A12) &\n" +
"                 - (r1215*ea5*d1rw*w4)/(r32*A32)                    &\n" +
"                 + (r15309*ea7*d1rw*w6)/(r128*A52)\n" +
"        !\n" +
"        np2 = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)        &\n" +
"              + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)\n" +
"        !\n" +
"        !\n" +
"        d1rnp2 =   f12*(Nine*ea2*d1rw*w)         &\n" +
"                 - (r81*ea4*d1rw*w3)/(Four*A)    &\n" +
"                 + (r2187*ea6*d1rw*w5)/(r32*A2)  &\n" +
"                 - (r6561*ea8*d1rw*w7)/(r32*A3)\n" +
"        !\n" +
"        ! The first term is\n" +
"        !\n" +
"        t1    = f12*(np1*piexperf + np2*expei)\n" +
"        d1st1 = f12*(d1spiexperf*np1 + d1sexpei*np2)\n" +
"        d1rt1 = f12*(d1rnp2*expei + d1rpiexperf*np1 +  &\n" +
"                     d1rexpei*np2 + d1rnp1*piexperf)\n" +
"        !\n" +
"        ! The factors for the main polynomoal in w and their derivatives\n" +
"        !\n" +
"        f2    = (f12)*ea1*srpi*A / DHsbw12\n" +
"        f2d1  = - ea1*srpi*A / (Four*DHsbw32)\n" +
"        d1sf2 = d1sHsbw*f2d1\n" +
"        d1rf2 = d1rHsbw*f2d1\n" +
"        !\n" +
"        f3    = (f12)*ea2*A / DHsbw\n" +
"        f3d1  = - ea2*A / (Two*DHsbw2)\n" +
"        d1sf3 = d1sHsbw*f3d1\n" +
"        d1rf3 = d1rHsbw*f3d1\n" +
"        !\n" +
"        f4    =  ea3*srpi*(-f98 / Hsbw12     &\n" +
"                 + f14*A / DHsbw32)\n" +
"        f4d1  = ea3*srpi*((Nine/(Sixteen*Hsbw32))-   &\n" +
"                          (Three*A/(Eight*DHsbw52)))\n" +
"        d1sf4 = d1sHsbw*f4d1\n" +
"        d1rf4 = d1rHsbw*f4d1\n" +
"        !\n" +
"        f5    = ea4*(One/r128) * (-r144*(One/Hsbw)   &\n" +
"                 + r64*(One/DHsbw2)*A)\n" +
"        f5d1  = ea4*((f98/Hsbw2)-(A/DHsbw3))\n" +
"        d1sf5 = d1sHsbw*f5d1\n" +
"        d1rf5 = d1rHsbw*f5d1\n" +
"        !\n" +
"        f6    = ea5*(Three*srpi*(Three*DHsbw52*(Nine*Hsbw-Two*A) &\n" +
"                 + Four*Hsbw32*A2))                              &\n" +
"                 / (r32*DHsbw52*Hsbw32*A)\n" +
"        f6d1  = ea5*srpi*((r27/(r32*Hsbw52))-        &\n" +
"                    (r81/(r64*Hsbw32*A))-            &\n" +
"                    ((Fifteen*A)/(Sixteen*DHsbw72)))\n" +
"        d1sf6 = d1sHsbw*f6d1\n" +
"        d1rf6 = d1rHsbw*f6d1\n" +
"        !\n" +
"        f7    = ea6*(((r32*A)/DHsbw3                 &\n" +
"                 + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32\n" +
"        d1sf7 = ea6*(Three*(r27*d1sH*DHsbw4*Hsbw*s2 +           &\n" +
"                Eight*d1sHsbw*A*(Three*DHsbw4 - Four*Hsbw3*A) + &\n" +
"                r54*DHsbw4*s*(Hsbw - d1sHsbw*s)*H))/            &\n" +
"                (r32*DHsbw4*Hsbw3*A)\n" +
"        d1rf7 = ea6*d1rHsbw*((f94/Hsbw3)-((Three*A)/DHsbw4)     &\n" +
"                           -((r81*s2*H)/(Sixteen*Hsbw3*A)))\n" +
"        !\n" +
"        f8    = ea7*(-Three*srpi*(-r40*Hsbw52*A3                &\n" +
"                 +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2))) &\n" +
"                 / (r128 * DHsbw72*Hsbw52*A2)\n" +
"        f8d1  = ea7*srpi*((r135/(r64*Hsbw72)) + (r729/(r256*Hsbw32*A2))  &\n" +
"                         -(r243/(r128*Hsbw52*A))                         &\n" +
"                         -((r105*A)/(r32*DHsbw92)))\n" +
"        d1sf8 = d1sHsbw*f8d1\n" +
"        d1rf8 = d1rHsbw*f8d1\n" +
"        !\n" +
"        f9    = (r324*ea6*eb1*DHsbw4*Hsbw*A                      &\n" +
"                + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2       &\n" +
"                + r324*Hsbw*A - r288*A2))) / (r128*DHsbw4*Hsbw3*A2)\n" +
"        f9d1  = -((r81*ea6*eb1)/(Sixteen*Hsbw3*A))               &\n" +
"                + ea8*((r27/(Four*Hsbw4))+(r729/(r128*Hsbw2*A2)) &\n" +
"                      -(r81/(Sixteen*Hsbw3*A))                   &\n" +
"                      -((r12*A/DHsbw5)))\n" +
"        d1sf9 = d1sHsbw*f9d1\n" +
"        d1rf9 = d1rHsbw*f9d1\n" +
"        !\n" +
"        t2t9    = f2*w  + f3*w2 + f4*w3 + f5*w4 + f6*w5          &\n" +
"                        + f7*w6 + f8*w7 + f9*w8\n" +
"        d1st2t9 = d1sf2*w + d1sf3*w2 + d1sf4*w3 + d1sf5*w4       &\n" +
"                          + d1sf6*w5 + d1sf7*w6 + d1sf8*w7       &\n" +
"                          + d1sf9*w8\n" +
"        d1rt2t9 = d1rw*f2 + d1rf2*w + Two*d1rw*f3*w   &\n" +
"                  + d1rf3*w2 + Three*d1rw*f4*w2       &\n" +
"                  + d1rf4*w3 + Four*d1rw*f5*w3        &\n" +
"                  + d1rf5*w4 + Five*d1rw*f6*w4        &\n" +
"                  + d1rf6*w5 + Six*d1rw*f7*w5         &\n" +
"                  + d1rf7*w6 + Seven*d1rw*f8*w6       &\n" +
"                  + d1rf8*w7 + Eight*d1rw*f9*w7 + d1rf9*w8\n" +
"        !\n" +
"        ! The final value of term1 for 0 < omega < wcutoff is:\n" +
"        !\n" +
"        term1 = t1 + t2t9 + t10\n" +
"        !\n" +
"        d1sterm1 = d1st1 + d1st2t9 + d1st10\n" +
"        d1rterm1 = d1rt1 + d1rt2t9 + d1rt10\n" +
"        !\n" +
"        ! The final value for the enhancement factor and its\n" +
"        ! derivatives is:\n" +
"        !\n" +
"        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)\n" +
"        !\n" +
"        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3    &\n";
      String code2="                              + d1sterm4 + d1sterm5)\n" +
"        !\n" +
"        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)\n" +
"        !\n" +
"      ENDIF\n" +
"\n" +
"END SUBROUTINE wpbe_analy_erfc_approx_grad\n" +
"!\n" +
"!---------------------------------------------------------------------\n" +
"function qe_erf(x)                      !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------------\n" +
"  !     Error function - computed from the rational approximations of\n" +
"  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.\n" +
"  !\n" +
"  !     for abs(x) le 0.47 erf is calculated directly\n" +
"  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)\n" +
"  USE kinds,   ONLY: DP\n" +
"  implicit none\n" +
"  REAL(DP), intent(in) :: x\n" +
"  REAL(DP) :: x2, p1 (4), q1 (4)\n" +
"  REAL(DP) :: qe_erf    !<GPU:qe_erf=>qe_erf_d>\n" +
"  data p1 / 2.426679552305318E2, 2.197926161829415E1, &\n" +
"            6.996383488619136d0,  -3.560984370181538E-2 /\n" +
"  data q1 / 2.150588758698612E2, 9.116490540451490E1, &\n" +
"            1.508279763040779E1, 1.000000000000000d0 /\n" +
"  !\n" +
"  if (abs (x) > 6.0d0) then\n" +
"     !\n" +
"     !  erf(6)=1-10^(-17) cannot be distinguished from 1\n" +
"     !\n" +
"     qe_erf = sign (1.0d0, x)                                                  !<GPU:qe_erf=>qe_erf_d>\n" +
"  else\n" +
"     if (abs (x)  <= 0.47d0) then\n" +
"        x2 = x**2\n" +
"        qe_erf=x *(p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &  !<GPU:qe_erf=>qe_erf_d>\n" +
"                / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )\n" +
"     else\n" +
"        qe_erf = 1.0d0 - qe_erfc(x)                                            !<GPU:qe_erf=>qe_erf_d,qe_erfc=>qe_erfc_d>\n" +
"     endif\n" +
"  endif\n" +
"  !\n" +
"  return\n" +
"end function qe_erf\n" +
"!\n" +
"!---------------------------------------------------------------------\n" +
"function qe_erfc(x)                      !<GPU:DEVICE>\n" +
"  !---------------------------------------------------------------------\n" +
"  !\n" +
"  !     erfc(x) = 1-erf(x)  - See comments in erf\n" +
"  !\n" +
"  USE kinds,   ONLY: DP\n" +
"  implicit none\n" +
"  !\n" +
"  REAL(DP),intent(in) :: x\n" +
"  REAL(DP)            :: qe_erfc                                         !<GPU:qe_erfc=>qe_erfc_d>\n" +
"  REAL(DP) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1\n" +
"  !\n" +
"  data p2 / 3.004592610201616E2,  4.519189537118719E2, &\n" +
"            3.393208167343437E2,  1.529892850469404E2, &\n" +
"            4.316222722205674E1,  7.211758250883094d0,   &\n" +
"            5.641955174789740E-1,-1.368648573827167E-7 /\n" +
"  data q2 / 3.004592609569833E2,  7.909509253278980E2, &\n" +
"            9.313540948506096E2,  6.389802644656312E2, &\n" +
"            2.775854447439876E2,  7.700015293522947E1, &\n" +
"            1.278272731962942E1,  1.000000000000000d0 /\n" +
"  data p3 /-2.996107077035422E-3,-4.947309106232507E-2, &\n" +
"           -2.269565935396869E-1,-2.786613086096478E-1, &\n" +
"           -2.231924597341847E-2 /\n" +
"  data q3 / 1.062092305284679E-2, 1.913089261078298E-1, &\n" +
"            1.051675107067932d0,    1.987332018171353d0,    &\n" +
"            1.000000000000000d0 /\n" +
"\n" +
"  data pim1 / 0.56418958354775629d0 /\n" +
"  !        ( pim1= sqrt(1/pi) )\n" +
"  ax = abs (x)\n" +
"  if (ax > 26.0d0) then\n" +
"     !\n" +
"     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);\n" +
"     !\n" +
"     qe_erfc = 0.0d0                                                            !<GPU:qe_erfc=>qe_erfc_d>\n" +
"  elseif (ax > 4.0d0) then\n" +
"     x2 = x**2\n" +
"     xm2 = (1.0d0 / ax) **2\n" +
"     qe_erfc = (1.0d0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &             !<GPU:qe_erfc=>qe_erfc_d>\n" +
"          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &\n" +
"          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &\n" +
"          (q3 (4) + xm2 * q3 (5) ) ) ) ) )\n" +
"  elseif (ax > 0.47d0) then\n" +
"     x2 = x**2\n" +
"     qe_erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &            !<GPU:qe_erfc=>qe_erfc_d>\n" +
"          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &\n" +
"          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &\n" +
"          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &\n" +
"          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )\n" +
"  else\n" +
"     qe_erfc = 1.0d0 - qe_erf(ax)                          !<GPU:qe_erfc=>qe_erfc_d, qe_erf=>qe_erf_d>\n" +
"  endif\n" +
"  !\n" +
"  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)\n" +
"  !\n" +
"  if (x < 0.0d0) qe_erfc = 2.0d0 - qe_erfc                                      !<GPU:qe_erfc=>qe_erfc_d>\n" +
"  !\n" +
"  return\n" +
"end function qe_erfc\n" +
"!\n" +
"\n" +
"!------------------------------------------------------------------\n" +
"FUNCTION EXPINT(n, x)                     !<GPU:DEVICE>\n" +
"!-----------------------------------------------------------------------\n" +
"! Evaluates the exponential integral E_n(x)\n" +
"! Parameters: maxit is the maximum allowed number of iterations,\n" +
"! eps is the desired relative error, not smaller than the machine precision,\n" +
"! big is a number near the largest representable floating-point number,\n" +
"! Inspired from Numerical Recipes\n" +
"!\n" +
"      USE kinds,   ONLY: DP\n" +
"      IMPLICIT NONE\n" +
"      INTEGER, INTENT(IN) :: n\n" +
"      REAL(DP), INTENT(IN) :: x\n" +
"      REAL(DP) :: expint                             !<GPU:expint=>expint_d>\n" +
"      INTEGER, parameter :: maxit=200\n" +
"      REAL(DP), parameter :: eps=1E-12, big=huge(x)*eps\n" +
"      REAL(DP), parameter :: euler = 0.577215664901532860606512d0\n" +
"!     EPS=1E-9, FPMIN=1E-30\n" +
"\n" +
"      INTEGER :: i, nm1, k\n" +
"      REAL(DP) :: a,b,c,d,del,fact,h,iarsum\n" +
"\n" +
"      IF (.NOT. ((n >= 0).AND.(x >= 0.0).AND.((x > 0.0).OR.(n > 1)))) THEN\n" +
"         !CALL errore('expint','bad arguments', 1)\n" +
"         STOP\n" +
"      END IF\n" +
"\n" +
"      IF (n == 0) THEN\n" +
"         expint = exp(-x)/x                                             !<GPU:expint=>expint_d>\n" +
"         RETURN\n" +
"      END IF\n" +
"      nm1 = n-1\n" +
"      IF (x == 0.0d0) THEN\n" +
"         expint = 1.0d0/nm1                                             !<GPU:expint=>expint_d>\n" +
"      ELSE IF (x > 1.0d0) THEN\n" +
"         b = x+n\n" +
"         c = big\n" +
"         d = 1.0d0/b\n" +
"         h = d\n" +
"         DO i=1,maxit\n" +
"            a = -i*(nm1+i)\n" +
"            b = b+2.0d0\n" +
"            d = 1.0d0/(a*d+b)\n" +
"            c = b+a/c\n" +
"            del = c*d\n" +
"            h = h*del\n" +
"            IF (ABS(del-1.0d0) <= EPS) EXIT\n" +
"         END DO\n" +
"         IF (i > maxit) STOP !CALL errore('expint','continued fraction failed',1)\n" +
"         expint = h*EXP(-x)                                             !<GPU:expint=>expint_d>\n" +
"      ELSE\n" +
"         IF (nm1 /= 0) THEN\n" +
"            expint = 1.0d0/nm1                                          !<GPU:expint=>expint_d>\n" +
"         ELSE\n" +
"            expint = -LOG(x)-euler                                      !<GPU:expint=>expint_d>\n" +
"         END IF\n" +
"         fact = 1.0d0\n" +
"         do i=1,maxit\n" +
"            fact = -fact*x/i\n" +
"            IF (i /= nm1) THEN\n" +
"               del = -fact/(i-nm1)\n" +
"            ELSE\n" +
"\n" +
"               iarsum = 0.0d0\n" +
"               do k=1,nm1\n" +
"                  iarsum = iarsum + 1.0d0/k\n" +
"               end do\n" +
"\n" +
"               del = fact*(-LOG(x)-euler+iarsum)\n" +
"!               del = fact*(-LOG(x)-euler+sum(1.0d0/arth(1,1,nm1)))\n" +
"            END IF\n" +
"            expint = expint + del                                       !<GPU:expint=>expint_d>\n" +
"            IF (ABS(del) < ABS(expint)*eps) EXIT                        !<GPU:expint=>expint_d>\n" +
"         END DO\n" +
"         IF (i > maxit) STOP !CALL errore('expint','series failed',1)\n" +
"      END IF\n" +
"END FUNCTION EXPINT\n" +
"!\n" +
"END MODULE";
       myWriter.write(code1+code2);
      myWriter.close();
      System.out.println("Successfully wrote to the file.");
    } catch (IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }
  }
}