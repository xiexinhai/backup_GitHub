import sys
import os.path
sys.path.insert(0, '/home/xinhai/MagnetControl/adc')
import aio320ra
import time
import json

def read_aio(aio, channel_id):
    return aio.analog_read_volt(channel_id, aio.DataRate.DR_860SPS)

def read_aio_all(aio):
    ret = []
    for channel in range(32):
        ret.append(read_aio(aio,channel))
    return ret

if __name__ == "__main__":
    aio = aio320ra.AIO_32_0RA_IRC(0x49, 0x3e)
    for channel in range(32):
        print('CH{:d}: {:2.3f}V'.format(channel, read_aio(aio,channel)))

    for channel in range(16):
        print('CH{:d}-CH{:d}: {:2.3f}V'.format(channel+16,channel, read_aio(aio,channel+16)-read_aio(aio,channel)))
