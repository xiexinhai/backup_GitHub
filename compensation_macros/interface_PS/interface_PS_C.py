import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX


if __name__ == "__main__":
    #When rebot the server, the USB ID of the power supply may change!!!!!!
    with open("../USBID_order.json","r") as file_obj:
        USBID_order = json.load(file_obj)
        file_obj.close()

    comm = PMX()
    comm.open(USBID_order[2], '19200')
    comm.Clear()
    result, data = comm.CheckVersion()
    print(result, data)

    print('Welcome to the interface of the power supply!')
    print('Enter 1 to turn on the output.')
    print('Enter 2 to turn off the output.')
    print('Enter 3 to set the voltage value.')
    print('Enter 4 to set the current value.')
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
            volt = input('Please enter the voltage (0~20V): \n>')
            if float(volt) >= 0 and float(volt)<= 20:
                comm.SetVoltage(float(volt))
                print('Voltage output value set to be U='+str(float(volt))+'V')
            else:
                print('Voltage value error!')

        elif int(inNum) == 4:
            curr = input('Please enter the current (0~5A): \n>')
            if float(curr) >= 0 and float(curr)<= 5:
                comm.SetCurrent(float(curr))
                print('Current output value set to be I='+str(float(curr))+'A')
            else:
                print('Current value error!')
        else:
            print('Error! Please enter the correct integer number!')


    comm.Output("OFF")
    comm.close()



