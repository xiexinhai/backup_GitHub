import math

class B_left_right:
    def __init__(self):
        self.factor = 25. # general factor
        self.mu0 = 4*math.pi*math.pow(10,-3) # constant, unit in Gs*m/A
        self.I = 1.0 # current, unit in A
        self.x0 = 0.179 # coordinate of coil plane, unit in m
        self.h1 = 0.316 # length of the squared coil along y axis, unit in m
        self.h2 = 0.378 # length of the squared coil along z axis, unit in m
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

    def set_x0(self,x0):
        self.x0 = x0

    def set_h1(self,h1):
        self.h1 = h1

    def set_h2(self,h2):
        self.h2 = h2

    def print(self):
        print('x=',self.x)
        print('y=',self.y)
        print('z=',self.z)

    def magnetic_x(self): #magnetic of only one squared coil at x0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (h2/2-z) / ( (x-x0)**2+(z-h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) )
        part2 = (h1/2-y) / ( (x-x0)**2+(y-h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) )
        part3 = (h2/2+z) / ( (x-x0)**2+(z+h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        part4 = (h1/2+y) / ( (x-x0)**2+(y+h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2+part3+part4)
        return ret

    def magnetic_x_symm(self): #magnetic of only one squared coil at symmetric position -x0 
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = -self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (h2/2-z) / ( (x-x0)**2+(z-h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) )
        part2 = (h1/2-y) / ( (x-x0)**2+(y-h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) )
        part3 = (h2/2+z) / ( (x-x0)**2+(z+h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        part4 = (h1/2+y) / ( (x-x0)**2+(y+h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2+part3+part4)
        return ret


    def magnetic_y(self): #magnetic of only one squared coil at x0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (x-x0) / ( (x-x0)**2+(y-h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) )
        part2 = (x0-x) / ( (x-x0)**2+(y+h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2)
        return ret

    def magnetic_y_symm(self): #magnetic of only one squared coil at symmetric position -x0 
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = -self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (x-x0) / ( (x-x0)**2+(y-h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) )
        part2 = (x0-x) / ( (x-x0)**2+(y+h1/2)**2 ) * ( -(-(h2/2) + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) + (h2/2 + z)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2)
        return ret

    def magnetic_z(self): #magnetic of only one squared coil at x0
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (x-x0) / ( (x-x0)**2+(z-h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) )
        part2 = (x0-x) / ( (x-x0)**2+(z+h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2)
        return ret

    def magnetic_z_symm(self):  #magnetic of only one squared coil at symmetric position -x0 
        factor = self.factor
        mu0 = self.mu0
        I = self.I
        x0 = -self.x0
        h1 = self.h1
        h2 = self.h2
        x = self.x
        y = self.y
        z = self.z

        part1 = (x-x0) / ( (x-x0)**2+(z-h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + (-(h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + (-(h2/2) + z)**2) )
        part2 = (x0-x) / ( (x-x0)**2+(z+h2/2)**2 ) * ( -(-(h1/2) + y)/math.sqrt((x - x0)**2 + (-(h1/2) + y)**2 + ((h2/2) + z)**2) + (h1/2 + y)/math.sqrt((x - x0)**2 + ((h1/2) + y)**2 + ((h2/2) + z)**2) )
        ret = factor*I*mu0/math.pi/4.*(part1+part2)
        return ret

    def magnetic3D(self):
        ret = [self.magnetic_x()+self.magnetic_x_symm(),self.magnetic_y()+self.magnetic_y_symm(),self.magnetic_z()+self.magnetic_z_symm()]
        return ret











