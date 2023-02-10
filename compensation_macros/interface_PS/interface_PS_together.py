import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX
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

    print('Welcome to the interface of the power supply!\nThis macro is to control the three power supplies together!')
    print('Enter 1 to turn on the output.')
    print('Enter 2 to turn off the output.')
    print('Enter 3 to set the voltage values.')
    print('Enter 4 to set the current values.')
    print('Enter 0 to close the interface.')

    while True:
        inNum = input('Please enter the integer number (0~4): \n>')
        if int(inNum) == 0:
            print('Interface macro closed!')
            break

        elif int(inNum) == 1:
            comm.Output("ON")
            print('Power supply output ON!')

        elif int(inNum) == 2:
            print('Power supply output OFF!')
            comm.Output("OFF")

        elif int(inNum) == 3:
            Volt = input('Please enter the voltage (0~20V) splitted by comma ",". eg: 1.0,2.0,3.0: \n>')
            volt = [float(n) for n in Volt.split(",")]
            if volt[0] >= 0 and volt[0]<= 20 and volt[1] >= 0 and volt[1]<= 20 and volt[2] >= 0 and volt[2]<= 20:
                comm.SetVoltage(volt[0])
                comm2.SetVoltage(volt[1])
                comm3.SetVoltage(volt[2])
                print('Voltage output value set to be U1='+str(volt[0])+'V, U2='+str(volt[1])+'V, U3='+str(volt[2])+'V')
            else:
                print('Voltage value error!')

        elif int(inNum) == 4:
            Curr = input('Please enter the current (0~5A) splitted by comma ",". eg: 1.0,2.0,3.0: \n>')
            curr = [float(n) for n in Curr.split(",")]
            if curr[0] >= 0 and curr[0]<= 5 and curr[1] >= 0 and curr[1]<= 5 and curr[2] >= 0 and curr[2]<= 5:
                comm.SetCurrent(curr[0])
                comm.SetCurrent(curr[1])
                comm.SetCurrent(curr[2])
                print('Current output value set to be I1='+str(curr[0])+'A, I2='+str(curr[1])+'A, I3='+str(curr[2])+'A')
            else:
                print('Current value error!')
        else:
            print('Error! Please enter the correct integer number!')


    comm.Output("OFF")
    comm.close()

    comm1.Output("OFF")
    comm1.close()

    comm2.Output("OFF")
    comm2.close()

