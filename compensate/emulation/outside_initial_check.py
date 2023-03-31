import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')
sys.path.insert(0, '/home/xinhai/MagnetControl/adc')

#import aio320ra
#from PMX import PMX
import time
import datetime
import math
import json
import numpy as np

deltaT = 1

Bchx = 1 # LR coil
Bchy = 0 # UD coil
Bchz = 2 # FB coil

T1ch = 3
T2ch = 4

#definition of the matrix related to three coils, determined by the position of the outside sensor
const_Bout = np.mat([-20.,-25.,6.]).T

OutBx  = 0.853237108
OutBxy = -0.717304637
OutBxz = -0.010788633

OutByx = -0.242167464
OutBy  = 0.978787453
OutByz = -0.047495391

OutBzx = 0.007704043
OutBzy = 0.094762699
OutBz  = -0.381223172

matOut = 100*np.mat([[OutBx,OutByx,OutBzx],[OutBxy,OutBy,OutBzy],[OutBxz,OutByz,OutBz]])
invmatOut = np.linalg.inv(matOut)

#definition of the matrix related to three coils, determined by the position of the center (the conpensate position)
const_Bcen = np.mat([-19.,-23.,5.]).T

CenBx  = 0.6335030219254629
CenBxy = 0
CenBxz = 0

CenByx = 0
CenBy  = 0.6585812023055216
CenByz = 0

CenBzx = 0
CenBzy = 0
CenBz  = -0.6530124174708793

matCen = 100*np.mat([[CenBx,CenByx,CenBzx],[CenBxy,CenBy,CenBzy],[CenBxz,CenByz,CenBz]])
invmatCen = np.linalg.inv(matCen)


def read_aio(aio, channel_id):
    return aio.analog_read_volt(channel_id, aio.DataRate.DR_860SPS)

def read_aio_all(aio):
    ret = []
    for channel in range(32):
        ret.append(read_aio(aio,channel))
    return ret

def check_safety(deltaBx,deltaBy,deltaBz,T1,T2,Ix,Iy,Iz):
    ret = True
    print('------Safety check------')
    print('Temperature:')
    print('T1 = ',T1,"degree")
    print('T2 = ',T2,"degree")

    print('Expected compensation B field at the center:')
    print('deltaBx = ',deltaBx,"muT")
    print('deltaBy = ',deltaBy,"muT")
    print('deltaBz = ',deltaBz,"muT")

    print('Needed current:')
    print('deltaIx(I_LR) = ',Ix,"A")
    print('deltaIy(I_UD) = ',Iy,"A")
    print('deltaIz(I_FB) = ',Iz,"A")

    #Current check
    if Ix > 5. or Iy > 5. or Iz > 5.:
        ret = False
        print("Needed current too large! Safety check failed.")


    #Temperature check
    if T1 > 30. or T2 > 30.:
        ret = False
        print("Temperature too high! Safety check failed.")

    #Probe signal check 
    Bx = abs(deltaBx)
    By = abs(deltaBy)
    Bz = abs(deltaBz)
    if Bx > 100 or By > 100 or Bz > 100:
        ret = False
        print("Detected magnetic field too large! Safety check failed.")

    if ret == True:
        print('Safety check is passed.')
    print('------------------------')
    return ret

def get_T(ADCreadout):
    T1 = ADCreadout[T1ch]
    T2 = ADCreadout[T2ch]
    return (T1,T2)

def get_sensorB(ADCreadout):
    #need to be modified when upgrading probes!!!
    #Vx -> By, Vy -> -Bx, Vz -> Bz

    Fx = -1*35*(ADCreadout[Bchx])
    Fy = 35*(ADCreadout[Bchy])
    Fz = 35*(ADCreadout[Bchz])

    return (Fx, Fy, Fz)


def get_deltaB(Fx, Fy, Fz, Ixnow, Iynow, Iznow):
    #need to be modified when upgrading probes!!!

    Inow = np.mat([Ixnow, Iynow, Iznow]).T
    F0out = np.mat([Fx, Fy, Fz]).T - np.dot(matOut,Inow)

    Fcen = np.dot(matCen,Inow) + F0out + const_Bcen - const_Bout

    deltaBx = Fcen[0,0]
    deltaBy = Fcen[1,0]
    deltaBz = Fcen[2,0]

    return (deltaBx, deltaBy, deltaBz)

def cal_deltaI(deltaBx,deltaBy,deltaBz):
    #need to be modified when selecting different compensation coordinate!!!
    #initial: Bx0 ~ , By0 ~ , Bz0 ~ 

    deltaB = -1*np.mat([deltaBx,deltaBy,deltaBz]).T

    deltaI = np.dot(invmatCen,deltaB)

    deltaIx = deltaI[0,0]
    deltaIy = deltaI[1,0]
    deltaIz = deltaI[2,0]

    return (deltaIx,deltaIy,deltaIz)

if __name__ == "__main__":
    print('Ambient field:')
    print('Bout = ('+str(const_Bout[0,0])+","+str(const_Bout[1,0])+","+str(const_Bout[2,0])+') muT')
    print('Bcenter = ('+str(const_Bcen[0,0])+","+str(const_Bcen[1,0])+","+str(const_Bcen[2,0])+') muT')

    ADCreadout = [0]*32
    ADCreadout[0] = const_Bout[1,0]/35.
    ADCreadout[1] = -1*const_Bout[0,0]/35.
    ADCreadout[2] = const_Bout[2,0]/35.
    ADCreadout[3] = 20.
    ADCreadout[4] = 20.

    Ixnow = 0.0 
    Iynow = 0.0
    Iznow = 0.0

    Fx, Fy, Fz = get_sensorB(ADCreadout)
    T1, T2 = get_T(ADCreadout)

    deltaBx, deltaBy, deltaBz = get_deltaB(Fx, Fy, Fz, Ixnow, Iynow, Iznow)

    deltaIx, deltaIy, deltaIz = cal_deltaI(deltaBx,deltaBy,deltaBz)

    Ix = Ixnow + deltaIx
    Iy = Iynow + deltaIy
    Iz = Iznow + deltaIz

    check_safety(deltaBx,deltaBy,deltaBz,T1,T2,Ix,Iy,Iz)

    print("--------------------Ctrl+C to stop------------------------")
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))
    print("Bx_sensor =",Fx,"muT")
    print("By_sensor =",Fy,"muT")
    print("Bz_sensor =",Fz,"muT")
    print("------------------------------")
    print("Ixnow =",Ixnow,"A")
    print("Iynow =",Iynow,"A")
    print("Iznow =",Iznow,"A")
    print("------------------------------")
    print("deltaIx =",deltaIx,"A")
    print("deltaIy =",deltaIy,"A")
    print("deltaIz =",deltaIz,"A")
    print("------------------------------")
    print("Ix =",Ix,"A")
    print("Iy =",Iy,"A")
    print("Iz =",Iz,"A")
    print("----------------------------------------------------------")
    Bxcen_after = const_Bcen + np.dot(matCen,np.mat([Ix,Iy,Iz]).T)
    print("After apply Ix Iy Iz, it should be:")
    print("Center field:")
    print("Bxcen_after =",Bxcen_after[0,0],"muT")
    print("Bycen_after =",Bxcen_after[1,0],"muT")
    print("Bzcen_after =",Bxcen_after[2,0],"muT")

    Bxout_after = const_Bout + np.dot(matOut,np.mat([Ix,Iy,Iz]).T)
    print("After apply Ix Iy Iz, it should be:")
    print("Center field:")
    print("Bxcen_after =",Bxout_after[0,0],"muT")
    print("Bycen_after =",Bxout_after[1,0],"muT")
    print("Bzcen_after =",Bxout_after[2,0],"muT")
















