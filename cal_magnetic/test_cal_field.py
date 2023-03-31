from B_left_right import B_left_right
from B_forward_back import B_forward_back
from B_up_down import B_up_down
import math

#def Bz(I=1): #magnetic field at center
#    factor = 29.955
#    mu0 = 4*math.pi*math.pow(10,-3) # constant, unit in Gs*m/A
#    z0 = 0.209
#    r0 = 0.22
#    num = I*mu0*r0*r0
#    den = math.pow(z0**2 + r0**2,3/2)
#    ret = factor*num/den
#    return ret
#
#def By(I=1): #magnetic field at center
#    factor = 25.077
#    mu0 = 4*math.pi*math.pow(10,-3) # constant, unit in Gs*m/A
#    y0 = 0.179
#    r0 = 0.209
#    num = I*mu0*r0*r0
#    den = math.pow(y0**2 + r0**2,3/2)
#    ret = factor*num/den
#    return ret
#
#def Bx(I=1): #magnetic field at center
#    factor = 25.256
#    mu0 = 4*math.pi*math.pow(10,-3) # constant, unit in Gs*m/A
#    x0 = 0.179
#    h1 = 0.316
#    h2 = 0.373
#    part1 = I*mu0*h1*h2/2/math.pi/(h2**2/4 + x0**2)/math.sqrt(h1**2/4 + h2**2/4 + x0**2)
#    part2 = I*mu0*h1*h2/2/math.pi/(h1**2/4 + x0**2)/math.sqrt(h1**2/4 + h2**2/4 + x0**2)
#
#    ret = factor*(part1+part2)
#    return ret


if __name__ == "__main__":
#    print("Center magnetic field when I = 1A:")
#    print("Bx=",Bx())
#    print("By=",By())
#    print("Bz=",Bz())
#
#    a = B_left_right()
#    a.set_I(1)
#    a.set_xyz(0.0,0.0,0.0)
#    print(a.magnetic3D())
#
#    b = B_up_down()
#    b.set_I(1)
#    b.set_xyz(0.0,0.0,0.0)
#    print(b.magnetic3D())
#
#    c = B_forward_back()
#    c.set_I(1)
#    c.set_xyz(0.0,0.0,0.0)
#    print(c.magnetic3D())


#    zlist1 = [-5.5, -20.5, -40.5, -90.5, -140.5, -180.5, -220.5]
#    zlist2 = [-33.5 ,-48.5 ,-68.5 ,-118.5 ,-168.5 ,-208.5 ,-248.5 ]
#    zlist3 = [-19.5 ,-34.5 ,-54.5 ,-104.5 ,-154.5 ,-194.5 ,-234.5 ]


#    a = B_up_down()
#    a.set_I(1)
#    for i in range(7):
#        #a.set_xyz(0.0,0.014,zlist1[i]/1000)
#        #print(a.magnetic3D()[0]*100)
#
#        a.set_xyz(-0.014,0.0,zlist2[i]/1000)
#        print(a.magnetic3D()[1]*100)
#
#        #a.set_xyz(0.014,-0.014,zlist3[i]/1000)
#        #print(a.magnetic3D()[2]*100)


#    a = B_left_right()
#    a.set_I(1)
#    for i in range(7):
#        a.set_xyz(0.0,0.014,zlist1[i]/1000)
#        print(a.magnetic3D()[0]*100)
#
#        #a.set_xyz(-0.014,0.0,zlist2[i]/1000)
#        #print(a.magnetic3D()[1]*100)
#
#        #a.set_xyz(0.014,-0.014,zlist3[i]/1000)
#        #print(a.magnetic3D()[2]*100)


#    a = B_forward_back()
#    a.set_I(1)
#    for i in range(7):
#        #a.set_xyz(0.0,0.014,zlist1[i]/1000)
#        #print(a.magnetic3D()[0]*100)
#
#        #a.set_xyz(-0.014,0.0,zlist2[i]/1000)
#        #print(a.magnetic3D()[1]*100)
#
#        a.set_xyz(0.014,-0.014,zlist3[i]/1000)
#        print(a.magnetic3D()[2]*100)

#    zlist1 = [0.207,0.187,0.167,0.147,0.127,0.107,0.087,0.067]
#    zlist2 = [0.157,0.137,0.117,0.097,0.077,0.057,0.037,0.017]
#    zlist3 = [0.182,0.162,0.142,0.122,0.102,0.082,0.062,0.042]

#    zlist1 = [0.2428,0.2228,0.2028,0.1828,0.1628,0.1428,0.1228,0.1028]
#    zlist2 = [0.2018,0.1818,0.1618,0.1418,0.1218,0.1018,0.0818,0.0618]
#    zlist3 = [0.223,0.203,0.183,0.163,0.143,0.123,0.103,0.083]

    #a = B_up_down()
    #a = B_left_right()
#    a = B_forward_back()
#    a.set_I(1)
#    for i in range(8):
#        #a.set_xyz(0.0,0.0,zlist1[i])
#        #print(a.magnetic3D()[0]*100)
#
#        #a.set_xyz(0.0,0.0,zlist2[i])
#        #print(a.magnetic3D()[1]*100)
#
#        a.set_xyz(0.0,0.0,zlist3[i])
#        print(a.magnetic3D()[2]*100)

    #a = B_up_down() #y-axis
    #a = B_left_right() #x-axis
    #a = B_forward_back() #z-axis

    #old
    #a.set_xyz(0,-0.014,-0.2585-0.014) #FM-3500: x-probe
    #a.set_xyz(-0.014,0,-0.2585+0.014) #FM-3500: y-probe
    #a.set_xyz(0.014,0.014,-0.2585) #FM-3500: z-probe

    #new
    #a.set_xyz(0,-0.014,-0.2745-0.014) #FM-3500: x-probe
    #a.set_xyz(-0.014,0,-0.2745+0.014) #FM-3500: y-probe
    #a.set_xyz(0.014,0.014,-0.2745) #FM-3500: z-probe


    #print(a.magnetic3D())

    #old
    #a.set_xyz(0,0,0.067) #FLC3-70: x-probe
    #a.set_xyz(0,0,0.017) #FLC3-70: y-probe
    #a.set_xyz(0,0,0.042) #FLC3-70: z-probe

    #new
    #a.set_xyz(0,0,0.1028) #FLC3-70: x-probe
    #a.set_xyz(0,0,0.0618) #FLC3-70: y-probe
    #a.set_xyz(0,0,0.083) #FLC3-70: z-probe

    #print(a.magnetic3D())

    #a = B_left_right() #x-axis
    #a = B_up_down() #y-axis
    a = B_forward_back() #z-axis

    #02.02 FM-3500 at center
    #a.set_xyz(0,0.014,-0.0315)
    #a.set_xyz(0.014,0,-0.0035)
    #a.set_xyz(-0.014,-0.014,-0.0175)

    #new @ S2
    #a.set_xyz(0,-0.014,-0.024) #FM-3500: x-probe
    #a.set_xyz(0.014,0,-0.050) #FM-3500: y-probe
    #a.set_xyz(-0.014,0.014,-0.036) #FM-3500: z-probe

    #print(a.magnetic3D()) 

    #02.02 FLC3-70 at outside trace
    #a.set_xyz(-0.132,-0.132,-0.0038)
    #a.set_xyz(-0.132,-0.132,-0.0448)
    #a.set_xyz(-0.132,-0.132,-0.025)

    #S2 area: FLC3-70 at outside trace
    a.set_xyz(-0.132,-0.132,-0.0058)
    print(a.magnetic3D()) 
    a.set_xyz(-0.132,-0.132,-0.0468)
    print(a.magnetic3D()) 
    a.set_xyz(-0.132,-0.132,-0.027)
    print(a.magnetic3D()) 



    #a = B_left_right() #x-axis
    #b = B_up_down() #y-axis
    #c = B_forward_back() #z-axis


    #chamber center
    #a.set_xyz(0,0,0)
    #b.set_xyz(0,0,0)
    #c.set_xyz(0,0,0)

    #target centor
    #a.set_xyz(0,0,-0.032)
    #b.set_xyz(0,0,-0.032)
    #c.set_xyz(0,0,-0.032)

    #place FM-3500 at target centor
    #a.set_xyz(0,-0.014,-0.032-0.014)
    #b.set_xyz(-0.014,0,-0.032+0.014)
    #c.set_xyz(0.014,0.014,-0.032)

    #print(a.magnetic3D()) # 0.6263892915507522
    #print(b.magnetic3D()) # 0.6510311532127788
    #print(c.magnetic3D()) # -0.6678149417091388

#    #Environment 
#    bkgx = 27.04
#    bkgy = -23.81
#    bkgz = -5.66
#
#    a.set_I(-bkgx/0.6263892915507522/100)
#    b.set_I(-bkgy/0.6510311532127788/100)
#    c.set_I(-bkgz/-1/0.6678149417091388/100)
#
#    n = 50
#    for i in range(n):
#        pos = 500.*(i)/n
#        print(pos)
#    print('----')
#
#    for i in range(n):
#        pos = 500.*(i)/n
#        a.set_xyz(0,0,pos/1000.)
#        print(a.magnetic3D()[0]*100+bkgx)
#    print('----')
#
#    for i in range(n):
#        pos = 500.*(i)/n
#        b.set_xyz(0,0,pos/1000.)
#        print(b.magnetic3D()[1]*100+bkgy)
#    print('----')
#
#    for i in range(n):
#        pos = 500.*(i)/n
#        c.set_xyz(0,0,pos/1000.)
#        print(c.magnetic3D()[2]*100+bkgz)
#    print('----')
#
#    a.set_xyz(0,0,-0.032)
#    b.set_xyz(0,0,-0.032)
#    c.set_xyz(0,0,-0.032)
#
#    print(a.magnetic3D()) 
#    print(b.magnetic3D()) 
#    print(c.magnetic3D()) 

