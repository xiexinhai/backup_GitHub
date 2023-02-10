import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX
import time
import datetime
import math
import json

if __name__ == "__main__":
    #When rebot the server, the USB ID of the power supply may change!!!!!!
    with open("../USBID_order.json","r") as file_obj:
        USBID_order = json.load(file_obj)
        file_obj.close()

    comm = PMX()
    comm.open(USBID_order[0], '19200')
    comm.Clear()
    
    comm2 = PMX()
    comm2.open(USBID_order[1], '19200')
    comm2.Clear()

    comm3 = PMX()
    comm3.open(USBID_order[2], '19200')
    comm3.Clear()

    f = open('log_PS_status.csv', "a")
    stat = []

    stat.append(comm.ReadVoltage()[1])
    stat.append(comm.ReadCurrent()[1])
    stat.append(comm2.ReadVoltage()[1])
    stat.append(comm2.ReadCurrent()[1])
    stat.append(comm3.ReadVoltage()[1])
    stat.append(comm3.ReadCurrent()[1])

    stat.insert(0,datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))

    ret = str(stat).strip('[')
    ret = ret.strip(']')

    comm.close()
    comm2.close()
    comm3.close()



