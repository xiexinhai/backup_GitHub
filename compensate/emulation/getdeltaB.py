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
    ADCreadout[7]=2.502
    ADCreadout[8]=2.711
    ADCreadout[9]=-0.537

    print('ADC: ')
    print(ADCreadout)
    deltaBx, deltaBy, deltaBz = get_deltaB(ADCreadout)
    print('deltaB: ')
    print('deltaBx = ',deltaBx,"muT")
    print('deltaBy = ',deltaBy,"muT")
    print('deltaBz = ',deltaBz,"muT")




