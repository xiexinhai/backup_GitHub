import sys
sys.path.insert(0, '/home/shift/MagnetControl/fieldcontrol')

from PMX import PMX

import json

import argparse

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

    parser = argparse.ArgumentParser(description="simple scripts")
    parser.add_argument("-c", "--coil", default="none", dest="coil")
    parser.add_argument("-m", "--mode", default="none", dest="mode") #current or voltage or none
    parser.add_argument("-v", "--voltage", type=float, default=0., dest="voltage")
    parser.add_argument("-i", "--current", type=float, default=0., dest="current")
    parser.add_argument("-s", "--switch", default="none", dest="switch") #ON or OFF or none
    results = parser.parse_args()

    coil = results.coil
    mode = results.mode
    voltage = results.voltage
    current = results.current
    switch = results.switch

    if coil == "UD":
        print('Controling: Up-down coil')
        if mode == "current":
            print('Setting output current I=',current,'A')
            comm.SetCurrent(current)
        elif mode == "voltage":
            print('Setting output voltage U=',voltage,'V')
            comm.SetVoltage(voltage)
        else:
            print('Do not change output values.')

        if switch == "ON":
            print('Output turn ON!')
            comm.Output("ON")
        elif switch == "OFF":
            print('Output turn OFF!')
            comm.Output("OFF")
        else:
            print('Do not switch ON or OFF.')

    elif coil == "LR":
        print('Controling: Left-right coil')
        if mode == "current":
            print('Setting output current I=',current,'A')
            comm2.SetCurrent(current)
        elif mode == "voltage":
            print('Setting output voltage U=',voltage,'V')
            comm2.SetVoltage(voltage)
        else:
            print('Do not change output values.')

        if switch == "ON":
            print('Output turn ON!')
            comm2.Output("ON")
        elif switch == "OFF":
            print('Output turn OFF!')
            comm2.Output("OFF")
        else:
            print('Do not switch ON or OFF.')



    elif coil == "FB":
        print('Controling: Forward-backward coil')
        if mode == "current":
            print('Setting output current I=',current,'A')
            comm3.SetCurrent(current)
        elif mode == "voltage":
            print('Setting output voltage U=',voltage,'V')
            comm3.SetVoltage(voltage)
        else:
            print('Do not change output values.')

        if switch == "ON":
            print('Output turn ON!')
            comm3.Output("ON")
        elif switch == "OFF":
            print('Output turn OFF!')
            comm3.Output("OFF")
        else:
            print('Do not switch ON or OFF.')

    else:
        print('Please select coil name as UD, LR or FB!')


