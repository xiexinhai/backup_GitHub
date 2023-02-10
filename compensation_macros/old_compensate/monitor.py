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
    with open("ADCreadout.json", 'w') as file_obj:
        json.dump(read_aio_all(aio), file_obj)
        file_obj.close()

#    print(read_aio_all(aio))
#    try:
#        print("Begin field monitor!")
#        while True:
#            time.sleep(1)
#            with open("ADCreadout.json", 'w') as file_obj:
#                json.dump(read_aio_all(aio), file_obj)
#    except KeyboardInterrupt:
#        print("Stopped by KeyboardInterrupt!")
#    except Exception as e:
#        print(e)
