import sys
sys.path.insert(0, '/home/shift/MagnetControl/fieldcontrol')
sys.path.insert(0, '/home/shift/MagnetControl/adc')

import aio320ra
from PMX import PMX
import time
import datetime
import math
import json
import numpy as np


#definition for time interval, ADC channel, constant B-field
deltaT = 1

chVx = 0 
chVy = 1 
chVz = 2 

T1ch = 3
T2ch = 4


const_Bout = np.mat([0,0,0]).T
const_Bcen = np.mat([0,0,0]).T


#definition of the matrix related to three coils, determined by the position of the outside sensor
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


def get_sensorB(ADCreadout):
    #need to be modified when upgrading probes!!!
    #Vx -> By, Vy -> -Bx, Vz -> Bz

    Fx = -1*35*(ADCreadout[chVy])
    Fy = 35*(ADCreadout[chVx])
    Fz = 35*(ADCreadout[chVz])

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
    #When rebot the server, the USB ID of the power supply may change!!!!!!
    with open("../USBID_order.json","r") as file_obj:
        USBID_order = json.load(file_obj)
        file_obj.close()

    #UD coil
    comm = PMX()
    comm.open(USBID_order[0], '19200')
    comm.Clear()
    
    #LR coil
    comm2 = PMX()
    comm2.open(USBID_order[1], '19200')
    comm2.Clear()

    #FB coil
    comm3 = PMX()
    comm3.open(USBID_order[2], '19200')
    comm3.Clear()

    triggerPrint = False
    try:
        comm.SetCurrent(0.)
        comm.Output("ON")

        comm2.SetCurrent(0.)
        comm2.Output("ON")

        comm3.SetCurrent(0.)
        comm3.Output("ON")

        while True:

            #ADC readout, minimum readout for current used AIO-32/0RA-IRC is 0.000306259 V
            aio = aio320ra.AIO_32_0RA_IRC(0x49, 0x3e)
            ADCreadout = read_aio_all(aio)
        
            #current values for power supply
            Ixnow = comm.ReadCurrent()[1]
            Iynow = comm2.ReadCurrent()[1]
            Iznow = comm3.ReadCurrent()[1]

            Fx, Fy, Fz = get_sensorB(ADCreadout)

            deltaBx, deltaBy, deltaBz = get_deltaB(Fx, Fy, Fz, Ixnow, Iynow, Iznow)

            deltaIx, deltaIy, deltaIz = cal_deltaI(deltaBx,deltaBy,deltaBz)

            Ix = Ixnow + deltaIx
            Iy = Iynow + deltaIy
            Iz = Iznow + deltaIz

            if check_safety(deltaBx,deltaBy,deltaBz,Ix,Iy,Iz) == False:
                triggerPrint = True
                break

            comm.SetCurrent(Iy)
            comm2.SetCurrent(Ix)
            comm3.SetCurrent(Iz)

            print("--------------------Ctrl+C to stop------------------------")
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))
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

            #compensation frequency/time-interval!!!
            time.sleep(deltaT)


    except KeyboardInterrupt:
        print("Stopped by KeyboardInterrupt!")
        comm.Output("OFF")
        comm2.Output("OFF")
        comm3.Output("OFF")

        #comm.SetCurrent(0.)
        #comm2.SetCurrent(0.)
        #comm3.SetCurrent(0.)

    except Exception as e:
        print(e)
        comm.Output("OFF")
        comm2.Output("OFF")
        comm3.Output("OFF")

        #comm.SetCurrent(0.)
        #comm2.SetCurrent(0.)
        #comm3.SetCurrent(0.)


    if triggerPrint == True:
        print("Stopped by safety check!")

    comm.Output("OFF")
    comm2.Output("OFF")
    comm3.Output("OFF")

    #comm.SetCurrent(0.)
    #comm2.SetCurrent(0.)
    #comm3.SetCurrent(0.)

    comm.close()
    comm2.close()
    comm3.close()



