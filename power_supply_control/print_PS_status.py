import sys
sys.path.insert(0, '/home/shift/MagnetControl/fieldcontrol')

from PMX import PMX
import time
import datetime
import math
import json

if __name__ == "__main__":
    #When rebot the server, the USB ID of the power supply may change!!!!!!
    with open("/home/shift/PS_order/USBID_order.json","r") as file_obj:
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


    print('-----------------------------\n')
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f") + '\n')
    print('comm1 U: ')
    print(comm.ReadVoltage()[1])
    print('comm1 I: ')
    print(comm.ReadCurrent()[1])
    print("---" + '\n')

    print('comm2 U: ')
    print(comm2.ReadVoltage()[1])
    print('comm2 I: ')
    print(comm2.ReadCurrent()[1])
    print("---" + '\n')

    print('comm3 U: ')
    print(comm3.ReadVoltage()[1])
    print('comm3 I: ')
    print(comm3.ReadCurrent()[1])
    print("---" + '\n')

    comm.close()
    comm2.close()
    comm3.close()



