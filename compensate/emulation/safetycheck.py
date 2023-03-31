import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')
sys.path.insert(0, '/home/xinhai/MagnetControl/adc')

#import aio320ra
#from PMX import PMX
import time
import datetime
import math
import json

deltaT = 1

Bchx = 8 # LR coil
Bchy = 7 # UD coil
Bchz = 9 # FB coil

T1ch = 14
T2ch = 15

def read_aio(aio, channel_id):
    return aio.analog_read_volt(channel_id, aio.DataRate.DR_860SPS)

def read_aio_all(aio):
    ret = []
    for channel in range(32):
        ret.append(read_aio(aio,channel))
    return ret

def check_safety(ADCreadout,Ix,Iy,Iz):
    ret = True
    print('------Safety check------')
    print('Temperature:')
    print('T1 = ',ADCreadout[14],"degree")
    print('T2 = ',ADCreadout[15],"degree")

    print('Detected B field(magnitude):')
    print('Bx = ',10.*abs(ADCreadout[Bchx] - ADCreadout[Bchx+16]),"muT")
    print('By = ',10.*abs(ADCreadout[Bchy] - ADCreadout[Bchy+16]),"muT")
    print('Bz = ',10.*abs(ADCreadout[Bchz] - ADCreadout[Bchz+16]),"muT")

    print('Needed current:')
    print('deltaIx(I_LR) = ',Ix,"A")
    print('deltaIx(I_UD) = ',Iy,"A")
    print('deltaIx(I_FB) = ',Iz,"A")

    #Current check
    if Ix > 5. or Iy > 5. or Iz > 5.:
        ret = False
        print("Needed current too large! Safety check failed.")


    #Temperature check
    T1 = abs(ADCreadout[T1ch] - ADCreadout[T1ch+16])
    T2 = abs(ADCreadout[T2ch] - ADCreadout[T2ch+16])
    if T1 > 30. or T2 > 30.:
        ret = False
        print("Temperature too high! Safety check failed.")

    #Probe signal check (FM-3500)
    Bx = 10.*abs(ADCreadout[Bchx] - ADCreadout[Bchx+16])
    By = 10.*abs(ADCreadout[Bchy] - ADCreadout[Bchy+16])
    Bz = 10.*abs(ADCreadout[Bchz] - ADCreadout[Bchz+16])
    if Bx > 100 or By > 100 or Bz > 100:
        ret = False
        print("Detected magnetic field too large! Safety check failed.")

    if ret == True:
        print('Safety check is passed.')
    print('------------------------')
    return ret

def get_deltaB(ADCreadout):
    #need to be modified when upgrading probes!!!

    deltaBx = 10*(ADCreadout[Bchx] - ADCreadout[Bchx+16])
    deltaBy = -1*10*(ADCreadout[Bchy] - ADCreadout[Bchy+16])
    deltaBz = 10*(ADCreadout[Bchz] - ADCreadout[Bchz+16])

    return (deltaBx, deltaBy, deltaBz)

def cal_deltaI(deltaBx,deltaBy,deltaBz):
    #need to be modified when selecting different compensation coordinate!!!

    #magnetic field at the coordinate of censor FM-3500 placed at target center
    #x-probe (0,0.014,-0.0315)
    #y-probe (0.014,0,-0.0035)
    #z-probe (-0.014,-0.014,-0.0175)

    #initial: Bx0 ~ 27, By0 ~ -25, Bz0 ~ -5

    Bx  = -1 * 0.62432752
    Bxz = -1 * -0.003345779

    By  = 0.657061546
    Byz = -0.003588381

    Bz  = -1 * -0.65464895

    deltaIx = -deltaBx/Bx/100

    deltaIy = -deltaBy/By/100

    deltaIz = -deltaBz/Bz/100 - deltaIx*Bxz/Bz - deltaIy*Byz/Bz

    return (deltaIx,deltaIy,deltaIz)

if __name__ == "__main__":
    #ADC readout
    #aio = aio320ra.AIO_32_0RA_IRC(0x49, 0x3e)
    #ADCreadout = read_aio_all(aio)
    ADCreadout = [0]*32
    ADCreadout[7]=55.02 
    ADCreadout[8]=2.711
    ADCreadout[9]=-0.537
    ADCreadout[14]=18
    ADCreadout[15]=19

    print('ADC: ')
    print(ADCreadout)

    deltaBx, deltaBy, deltaBz = get_deltaB(ADCreadout)

    deltaIx, deltaIy, deltaIz = cal_deltaI(deltaBx,deltaBy,deltaBz)

    check_safety(ADCreadout,deltaIx, deltaIy, deltaIz)


