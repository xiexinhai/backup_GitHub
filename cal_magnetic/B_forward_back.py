import math
from scipy import integrate

class B_forward_back:
    def __init__(self):
        self.factor = -30 # general factor
        self.mu0 = 4*math.pi*math.pow(10,-3) # constant, unit in Gs*m/A
        self.I = 1.0 # current, unit in A
        self.z0 = 0.209 # coordinate of coil plane, unit in m
        self.r0 = 0.22 # length of the squared coil, unit in m
        #position unit of m
        self.x = 0.
        self.y = 0.
        self.z = 0.

    def set_xyz(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

    def set_factor(self,factor):
        self.factor = factor

    def set_mu0(self,mu0):
        self.mu0 = mu0

    def set_I(self,I):
        self.I = I

    def set_z0(self,z0):
        self.z0 = z0

    def set_r0(self,r0):
        self.r0 = r0

    def print(self):
        print('x=',self.x)
        print('y=',self.y)
        print('z=',self.z)

    def magnetic_x(self,th): #magnetic of coil pair at z0 and -z0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        z0 = self.z0
        r0 = self.r0
        x = self.x
        y = self.y
        z = self.z

        part1 = (z - z0)*math.cos(th)/math.pow((z - z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)
        part2 = (z + z0)*math.cos(th)/math.pow((z + z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)

        ret = factor*I*mu0*r0/4./math.pi*(part1+part2)
        return ret

    def magnetic_y(self,th): #magnetic of coil pair at z0 and -z0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        z0 = self.z0
        r0 = self.r0
        x = self.x
        y = self.y
        z = self.z

        part1 = (z - z0)*math.sin(th)/math.pow((z - z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)
        part2 = (z + z0)*math.sin(th)/math.pow((z + z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)

        ret = factor*I*mu0*r0/4./math.pi*(part1+part2)
        return ret

    def magnetic_z(self,th): #magnetic of coil pair at z0 and -z0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        z0 = self.z0
        r0 = self.r0
        x = self.x
        y = self.y
        z = self.z

        part1 = (r0 - x*math.cos(th) - y*math.sin(th))/math.pow((z - z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)
        part2 = (r0 - x*math.cos(th) - y*math.sin(th))/math.pow((z + z0)**2 + (x - r0*math.cos(th))**2 + (y - r0*math.sin(th))**2,1.5)

        ret = factor*I*mu0*r0/4./math.pi*(part1+part2)
        return ret

    def magnetic3D(self):
        retx,errx = integrate.quad(self.magnetic_x,0.,2*math.pi)
        rety,erry = integrate.quad(self.magnetic_y,0.,2*math.pi)
        retz,errz = integrate.quad(self.magnetic_z,0.,2*math.pi)
        ret = [retx,rety,retz]
        return ret











