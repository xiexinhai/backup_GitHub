import sys
import os.path
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX
import time


if __name__ == "__main__":
    comm = PMX()
    comm.open('/dev/ttyUSB0', '19200')
    comm.Clear()
    print("ID of comm1:")
    print(comm.CheckVersion()[1])

    comm2 = PMX()
    comm2.open('/dev/ttyUSB1', '19200')
    comm2.Clear()
    print("ID of comm2:")
    print(comm2.CheckVersion()[1])

    comm3 = PMX()
    comm3.open('/dev/ttyUSB2', '19200')
    comm3.Clear()
    print("ID of comm3:")
    print(comm3.CheckVersion()[1])

#    comm.SetVoltage(0.1)
#    comm.Output("ON")
#    time.sleep(1)
#    comm2.SetVoltage(0.1)
#    comm2.Output("ON")
#    time.sleep(1)
#    comm3.SetVoltage(0.1)
#    comm3.Output("ON")
#    time.sleep(1)

#    comm.SetCurrent(0.02)
#    comm.Output("ON")
#    time.sleep(1)
#    comm2.SetCurrent(0.02)
#    comm2.Output("ON")
#    time.sleep(1)
#    comm3.SetCurrent(0.02)
#    comm3.Output("ON")
#    time.sleep(1)

#    time.sleep(2)
#    comm2.SetCurrent(0.02)
#    comm2.Output("ON")
#    time.sleep(2)
#    comm2.SetCurrent(0.04)
#    comm2.Output("ON")
#    time.sleep(2)
#    comm2.SetCurrent(0.06)
#    comm2.Output("ON")
#    time.sleep(2)
#    comm2.Output("OFF")

    time.sleep(2)
    comm3.SetCurrent(0.02)
    comm3.Output("ON")
    time.sleep(2)
    comm3.SetCurrent(0.04)
    comm3.Output("ON")
    time.sleep(2)
    comm3.SetCurrent(0.06)
    comm3.Output("ON")
    time.sleep(2)
    comm3.Output("OFF")

    comm.close()
    comm2.close()
    comm3.close()



