"""
pmt_s1param.py
author: jrenner

Defines the PMT S1 parameterization functions as:

N(x; y,z) = sum(n)( sum(m)(c_nm*z^m) * y^n ) for x = 0 to x = 9

where (x,y,z) are the coordinates of the light emission point.

Because the response is characterized over several time bins, we
have several values for the coefficients.


Corona 0:

( ( 0.000847293862083*z^0 + 1.12820838693e-06*z^1 + 1.46367921997e-08*z^2 - 4.32332218662e-11*z^3 + 4.20390371592e-14*z^4)*y^0+ 
  (-4.13543734655e-06*z^0+7.32891262211e-08*z^1-4.66641479593e-10*z^2+1.28230905971e-12*z^3-1.21936525685e-15*z^4)*y^1+
  ( 2.24986625018e-08*z^0-3.72555799541e-10*z^1+2.3279640873e-12*z^2-6.4575693783e-15*z^3+6.49196050032e-18*z^4)*y^2)*(x==0)+

( ( 0.00084248854649*z^0+8.67182971828e-07*z^1+1.8338855268e-08*z^2-5.75433358781e-11*z^3+5.89000912166e-14*z^4)*y^0+
  (-4.33588677325e-06*z^0+8.69782164639e-08*z^1-6.09238978308e-10*z^2+1.78061797013e-12*z^3-1.78162711513e-15*z^4)*y^1+
  ( 2.29501934464e-08*z^0-4.36331390467e-10*z^1+3.03484558509e-12*z^2-8.95201123779e-15*z^3+9.26256015609e-18*z^4)*y^2)*(x==1)+

( ( 0.000732428873014*z^0+2.80919627426e-06*z^1+6.72731123532e-09*z^2-2.96443952299e-11*z^3+3.59847484462e-14*z^4)*y^0+
  (-1.6703799368e-06*z^0+3.9171805269e-08*z^1-3.14332697444e-10*z^2+1.03910136361e-12*z^3-1.1492098439e-15*z^4)*y^1+
  ( 8.5861961752e-09*z^0-1.77392729387e-10*z^1+1.38652801674e-12*z^2-4.59219535382e-15*z^3+5.25129433087e-18*z^4)*y^2)*(x==2)+
  
( ( 0.000752458415646*z^0+2.77252912512e-06*z^1+5.44509107657e-09*z^2-2.27925446997e-11*z^3+2.68066260142e-14*z^4)*y^0+
  (-2.33881442595e-06*z^0+4.5867008956e-08*z^1-3.29533569189e-10*z^2+9.99973743789e-13*z^3-1.05441230372e-15*z^4)*y^1+
  ( 1.17174974768e-08*z^0-2.12141628805e-10*z^1+1.47167723779e-12*z^2-4.33728144531e-15*z^3+4.54350739451e-18*z^4)*y^2)*(x==3)+

( ( 0.000776571854673*z^0+2.15010260373e-06*z^1+1.01649467099e-08*z^2-3.70271776251e-11*z^3+4.18687541532e-14*z^4)*y^0+
  (-1.34163957161e-06*z^0+2.62241832717e-08*z^1-1.98895873488e-10*z^2+6.53789435408e-13*z^3-7.65572090217e-16*z^4)*y^1+
  ( 3.3834426945e-09*z^0-4.33509789053e-11*z^1+3.16078045946e-13*z^2-1.13245162388e-15*z^3+1.5246164098e-18*z^4)*y^2)*(x==4)+
  
( ( 0.000767717377656*z^0+2.51845444939e-06*z^1+7.03364223155e-09*z^2-2.73835439434e-11*z^3+3.22797892474e-14*z^4)*y^0+
  (-1.220581228e-06*z^0+2.10093162961e-08*z^1-1.50856117174e-10*z^2+4.87558509345e-13*z^3-5.98577427426e-16*z^4)*y^1+
  ( 3.4091537289e-09*z^0-3.94535794693e-11*z^1+2.48603658047e-13*z^2-7.88789711394e-16*z^3+1.0741082416e-18*z^4)*y^2)*(x==5)+
  
( ( 0.000794960238901*z^0+2.12083284766e-06*z^1+8.9709860809e-09*z^2-3.16348278492e-11*z^3+3.60784640191e-14*z^4)*y^0+
  (-1.30835704374e-06*z^0+2.24186334326e-08*z^1-1.53497847752e-10*z^2+4.71019780999e-13*z^3-5.7892421648e-16*z^4)*y^1+
  ( 3.16706354563e-09*z^0-3.99036311342e-11*z^1+2.47103301703e-13*z^2-7.1407255215e-16*z^3+9.47783757899e-19*z^4)*y^2)*(x==6)+
  
( ( 0.000894021337629*z^0+2.93630824503e-07*z^1+2.01983730098e-08*z^2-5.91815916823e-11*z^3+6.00601865898e-14*z^4)*y^0+
  (-2.11097375466e-06*z^0+3.70695125768e-08*z^1-2.42586196621e-10*z^2+6.79638113206e-13*z^3-7.61607312097e-16*z^4)*y^1+
  ( 4.83200840806e-09*z^0-7.11478768834e-11*z^1+4.37737521491e-13*z^2-1.152562382e-15*z^3+1.31599177879e-18*z^4)*y^2)*(x==7)+
  
( ( 0.000847156563182*z^0+1.10634145132e-06*z^1+1.55872655232e-08*z^2-4.87932449839e-11*z^3+5.21144009566e-14*z^4)*y^0+
  (-1.72786872718e-06*z^0+3.03159307137e-08*z^1-2.01391242535e-10*z^2+5.72852531535e-13*z^3-6.68689768612e-16*z^4)*y^1+
  ( 4.20116643378e-09*z^0-6.0907291675e-11*z^1+3.74598963761e-13*z^2-9.75481047474e-16*z^3+1.14340564221e-18*z^4)*y^2)*(x==8)
  
Corona 1:

( ( 0.000841166133443*z^0-2.50662987502e-06*z^1+3.96165249784e-08*z^2-1.21941922328e-10*z^3+1.37622349169e-13*z^4)*y^0+
  (-2.01397956462e-05*z^0+4.45443850687e-07*z^1-3.0947334827e-09*z^2+8.0531377907e-12*z^3-6.97954517234e-15*z^4)*y^1+
  ( 6.63488920536e-07*z^0-1.47855311518e-08*z^1+1.03575846247e-10*z^2-2.74855446407e-13*z^3+2.43916666366e-16*z^4)*y^2)*(x==0)+

( ( 0.000671053377997*z^0+1.21056036512e-06*z^1+1.4028501197e-08*z^2-5.57154062485e-11*z^3+8.05197960169e-14*z^4)*y^0+
  (-3.54082414564e-06*z^0+3.44345823405e-08*z^1-1.7637222223e-10*z^2+6.40735775612e-13*z^3-8.75780985886e-16*z^4)*y^1+
  ( 3.06549018222e-07*z^0-4.77502837362e-09*z^1+2.95848614001e-11*z^2-8.55440232846e-14*z^3+8.9332526184e-17*z^4)*y^2)*(x==1)+

( ( 0.000582968586087*z^0+3.8243264999e-06*z^1-6.99971083498e-09*z^2+5.48977659768e-12*z^3+2.25448342599e-14*z^4)*y^0+
  ( 1.26560903267e-05*z^0-3.35410044923e-07*z^1+2.50153725071e-09*z^2-6.96654424608e-12*z^3+6.45333760581e-15*z^4)*y^1+
  (-2.40135921026e-07*z^0+6.36455779655e-09*z^1-4.65631983122e-11*z^2+1.27283214687e-13*z^3-1.1673427649e-16*z^4)*y^2)*(x==2)+

( ( 0.00056905431165*z^0+2.73235976055e-06*z^1+5.81350251173e-09*z^2-3.5453431775e-11*z^3+6.157115573e-14*z^4)*y^0+
  ( 7.82908474192e-06*z^0-1.39808673637e-07*z^1+8.67241120115e-10*z^2-2.26308472643e-12*z^3+2.10323452842e-15*z^4)*y^1+
  (-9.16513615287e-08*z^0+1.69329156649e-09*z^1-1.09103522349e-11*z^2+2.93812918693e-14*z^3-2.83328446528e-17*z^4)*y^2)*(x==3)+

( ( 0.000628184771147*z^0+1.8774773729e-06*z^1+1.02409451967e-08*z^2-4.56021904244e-11*z^3+7.04994088931e-14*z^4)*y^0+
  ( 1.18690495637e-06*z^0-2.30076470562e-08*z^1+1.46965359108e-10*z^2-3.8432647266e-13*z^3+3.40362658134e-16*z^4)*y^1+
  (-8.37054332346e-09*z^0+2.00172365974e-10*z^1-1.53099182625e-12*z^2+4.82966187083e-15*z^3-5.39602207694e-18*z^4)*y^2)*(x==4)+
  
( ( 0.000750079232543*z^0-7.75009879664e-07*z^1+2.88540978604e-08*z^2-9.75944608638e-11*z^3+1.21502868551e-13*z^4)*y^0+
  (-2.47293671567e-06*z^0+5.52087948981e-08*z^1-4.06163071363e-10*z^2+1.21812652847e-12*z^3-1.30122363773e-15*z^4)*y^1+
  ( 8.17249715316e-09*z^0-1.62656199423e-10*z^1+1.1141932498e-12*z^2-3.07723199194e-15*z^3+2.91883655113e-18*z^4)*y^2)*(x==5)+
  
( ( 0.000735109173329*z^0-4.52465761463e-07*z^1+2.68283918181e-08*z^2-9.32941797276e-11*z^3+1.20635848258e-13*z^4)*y^0+
  (-2.11869087455e-06*z^0+4.77074139908e-08*z^1-3.63760181212e-10*z^2+1.16098588878e-12*z^3-1.35467419855e-15*z^4)*y^1+
  ( 5.57100464613e-09*z^0-1.13342816042e-10*z^1+8.35840986049e-13*z^2-2.58443168804e-15*z^3+2.8567503534e-18*z^4)*y^2)*(x==6)+
  
( ( 0.000672270379336*z^0+1.02607902658e-06*z^1+1.61605482999e-08*z^2-6.22520636889e-11*z^3+9.07722806953e-14*z^4)*y^0+
  (-9.78957929879e-07*z^0+2.27936779553e-08*z^1-1.86939109714e-10*z^2+6.52238112102e-13*z^3-8.76558900574e-16*z^4)*y^1+
  ( 2.13536821414e-09*z^0-4.17040191573e-11*z^1+3.36915587974e-13*z^2-1.16675107176e-15*z^3+1.51633623612e-18*z^4)*y^2)*(x==7)+
  
( ( 0.000653141121787*z^0+1.58768761412e-06*z^1+1.15875295596e-08*z^2-4.67516073468e-11*z^3+7.42067525677e-14*z^4)*y^0+
  (-4.79910167246e-07*z^0+1.19890917766e-08*z^1-1.08077058814e-10*z^2+4.03522206326e-13*z^3-6.21161864206e-16*z^4)*y^1+
  ( 8.05608704706e-10*z^0-1.50360102487e-11*z^1+1.48437323618e-13*z^2-5.87986320847e-16*z^3+9.22272890447e-19*z^4)*y^2)*(x==8)  
  

"""
import numpy as np

#  
# Parameter table: indexed as
#
#  ptbl[corona][x][y][z]
#
# where corona runs from 0 to n_tbins-1
#            x runs from 0 to 9
#            y runs from 0 to 2
#            z runs from 0 to 4.
#
ptbl = [

        [
         [
          [ 0.000847293862083,  1.12820838693e-06,  1.46367921997e-08, -4.32332218662e-11,  4.20390371592e-14],
          [-4.13543734655e-06,  7.32891262211e-08, -4.66641479593e-10,  1.28230905971e-12, -1.21936525685e-15],
          [ 2.24986625018e-08, -3.72555799541e-10,  2.3279640873e-12,  -6.4575693783e-15,   6.49196050032e-18]
         ],
         [
          [ 0.00084248854649,   8.67182971828e-07,  1.8338855268e-08,  -5.75433358781e-11,  5.89000912166e-14],
          [-4.33588677325e-06,  8.69782164639e-08, -6.09238978308e-10,  1.78061797013e-12, -1.78162711513e-15],
          [ 2.29501934464e-08, -4.36331390467e-10,  3.03484558509e-12, -8.95201123779e-15,  9.26256015609e-18]
         ],
         [
          [ 0.000732428873014,  2.80919627426e-06,  6.72731123532e-09, -2.96443952299e-11,  3.59847484462e-14],
          [-1.6703799368e-06,   3.9171805269e-08,  -3.14332697444e-10,  1.03910136361e-12, -1.1492098439e-15],
          [ 8.5861961752e-09,  -1.77392729387e-10,  1.38652801674e-12, -4.59219535382e-15,  5.25129433087e-18]
         ],
         [
          [ 0.000752458415646,  2.77252912512e-06,  5.44509107657e-09, -2.27925446997e-11,  2.68066260142e-14],
          [-2.33881442595e-06,  4.5867008956e-08,  -3.29533569189e-10,  9.99973743789e-13, -1.05441230372e-15],
          [ 1.17174974768e-08, -2.12141628805e-10,  1.47167723779e-12, -4.33728144531e-15,  4.54350739451e-18]
         ],
         [
          [ 0.000776571854673,  2.15010260373e-06,  1.01649467099e-08, -3.70271776251e-11,  4.18687541532e-14],
          [-1.34163957161e-06,  2.62241832717e-08, -1.98895873488e-10,  6.53789435408e-13, -7.65572090217e-16],
          [ 3.3834426945e-09,  -4.33509789053e-11,  3.16078045946e-13, -1.13245162388e-15,  1.5246164098e-18]
         ],
         [
          [ 0.000767717377656,  2.51845444939e-06,  7.03364223155e-09, -2.73835439434e-11,  3.22797892474e-14],
          [-1.220581228e-06,    2.10093162961e-08, -1.50856117174e-10,  4.87558509345e-13, -5.98577427426e-16],
          [ 3.4091537289e-09,  -3.94535794693e-11,  2.48603658047e-13, -7.88789711394e-16,  1.0741082416e-18]
         ],
         [
          [ 0.000794960238901,  2.12083284766e-06,  8.9709860809e-09,  -3.16348278492e-11,  3.60784640191e-14],
          [-1.30835704374e-06,  2.24186334326e-08, -1.53497847752e-10,  4.71019780999e-13, -5.7892421648e-16],
          [ 3.16706354563e-09, -3.99036311342e-11,  2.47103301703e-13, -7.1407255215e-16,   9.47783757899e-19]
         ],
         [
          [ 0.000894021337629,  2.93630824503e-07,  2.01983730098e-08, -5.91815916823e-11,  6.00601865898e-14],
          [-2.11097375466e-06,  3.70695125768e-08, -2.42586196621e-10,  6.79638113206e-13, -7.61607312097e-16],
          [ 4.83200840806e-09, -7.11478768834e-11,  4.37737521491e-13, -1.152562382e-15,    1.31599177879e-18]
         ],
         [
          [ 0.000847156563182,  1.10634145132e-06,  1.55872655232e-08, -4.87932449839e-11,  5.21144009566e-14],
          [-1.72786872718e-06,  3.03159307137e-08, -2.01391242535e-10,  5.72852531535e-13, -6.68689768612e-16],
          [ 4.20116643378e-09, -6.0907291675e-11,   3.74598963761e-13, -9.75481047474e-16,  1.14340564221e-18]
         ]
        ],
        
        
        [
         [
          [ 0.000841166133443, -2.50662987502e-06,  3.96165249784e-08, -1.21941922328e-10,  1.37622349169e-13],
          [-2.01397956462e-05,  4.45443850687e-07, -3.0947334827e-09,   8.0531377907e-12,  -6.97954517234e-15],
          [ 6.63488920536e-07, -1.47855311518e-08,  1.03575846247e-10, -2.74855446407e-13,  2.43916666366e-16]
         ],
         [
          [ 0.000671053377997,  1.21056036512e-06,  1.4028501197e-08,  -5.57154062485e-11,  8.05197960169e-14],
          [-3.54082414564e-06,  3.44345823405e-08, -1.7637222223e-10,   6.40735775612e-13, -8.75780985886e-16],
          [ 3.06549018222e-07, -4.77502837362e-09,  2.95848614001e-11, -8.55440232846e-14,  8.9332526184e-17]
         ],
         [
          [ 0.000582968586087,  3.8243264999e-06,  -6.99971083498e-09,  5.48977659768e-12,  2.25448342599e-14],
          [ 1.26560903267e-05, -3.35410044923e-07,  2.50153725071e-09, -6.96654424608e-12,  6.45333760581e-15],
          [-2.40135921026e-07,  6.36455779655e-09, -4.65631983122e-11,  1.27283214687e-13, -1.1673427649e-16]
         ],
         [
          [ 0.00056905431165,   2.73235976055e-06,  5.81350251173e-09, -3.5453431775e-11,   6.157115573e-14],
          [ 7.82908474192e-06, -1.39808673637e-07,  8.67241120115e-10, -2.26308472643e-12,  2.10323452842e-15],
          [-9.16513615287e-08,  1.69329156649e-09, -1.09103522349e-11,  2.93812918693e-14, -2.83328446528e-17]
         ],
         [
          [ 0.000628184771147,  1.8774773729e-06,   1.02409451967e-08, -4.56021904244e-11,  7.04994088931e-14],
          [ 1.18690495637e-06, -2.30076470562e-08,  1.46965359108e-10, -3.8432647266e-13,   3.40362658134e-16],
          [-8.37054332346e-09,  2.00172365974e-10, -1.53099182625e-12,  4.82966187083e-15, -5.39602207694e-18]
         ],
         [
          [ 0.000750079232543, -7.75009879664e-07,  2.88540978604e-08, -9.75944608638e-11,  1.21502868551e-13],
          [-2.47293671567e-06,  5.52087948981e-08, -4.06163071363e-10,  1.21812652847e-12, -1.30122363773e-15],
          [ 8.17249715316e-09, -1.62656199423e-10,  1.1141932498e-12,  -3.07723199194e-15,  2.91883655113e-18]
         ],
         [
          [ 0.000735109173329, -4.52465761463e-07,  2.68283918181e-08, -9.32941797276e-11,  1.20635848258e-13],
          [-2.11869087455e-06,  4.77074139908e-08, -3.63760181212e-10,  1.16098588878e-12, -1.35467419855e-15],
          [ 5.57100464613e-09, -1.13342816042e-10,  8.35840986049e-13, -2.58443168804e-15,  2.8567503534e-18]
         ],
         [
          [ 0.000672270379336,  1.02607902658e-06,  1.61605482999e-08, -6.22520636889e-11,  9.07722806953e-14],
          [-9.78957929879e-07,  2.27936779553e-08, -1.86939109714e-10,  6.52238112102e-13, -8.76558900574e-16],
          [ 2.13536821414e-09, -4.17040191573e-11,  3.36915587974e-13, -1.16675107176e-15,  1.51633623612e-18]
         ],
         [
          [ 0.000653141121787,  1.58768761412e-06,  1.15875295596e-08, -4.67516073468e-11,  7.42067525677e-14],
          [-4.79910167246e-07,  1.19890917766e-08, -1.08077058814e-10,  4.03522206326e-13, -6.21161864206e-16],
          [ 8.05608704706e-10, -1.50360102487e-11,  1.48437323618e-13, -5.87986320847e-16,  9.22272890447e-19]
         ]
        ]
        
       ]
       
# Dimensions of ptbl
n_tbins = len(ptbl)
n_xdim = len(ptbl[0])
n_ydim = len(ptbl[0][0])
n_zdim = len(ptbl[0][0][0])
#print "Dimensions (tbins,x,y,z) = ({0},{1},{2},{3})".format(n_tbins,n_xdim,n_ydim,n_zdim)

# Maximum extent of dimensions
xmin = 0.; xmax = 8
ymin = 0.; ymax = 217.5
zmin = 65.; zmax = 475.

# Return the SiPM response for the specified time bin and (x,y,z).
# Note that x must be an integer from 0 to 8, and y and z are continuous.
def pmt_s1par(tbin,x,y,z):

    # Ensure the input values are valid
    if(tbin < 0 or tbin >= n_tbins):
        print "Invalid time bin in sipm_param: returning 0.0 ..."
        return 0.0
    if(not isinstance(x,int) or x < xmin or x > xmax):
        print "Invalid x value: must be integer in range [{0},{1}]".format(xmin,xmax)
        return 0.0
    if(y < ymin or y > ymax):
        print "Invalid y value: must be in range [{0},{1}]".format(ymin,ymax)
        return 0.0
    if(z < zmin or z > zmax):
        print "Invalid z value: must be in range [{0},{1}]".format(zmin,zmax)
        return 0.0
        
    # Calculate the response based on the parametrization.
    vpar = 0
    for yi in range(n_ydim):
            
        zpar = 0
        for zi in range(n_zdim):
            zpar += ptbl[tbin][x][yi][zi]*z**zi
            
        vpar += zpar*y**yi
        
            
    return vpar