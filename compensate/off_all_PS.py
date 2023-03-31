import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX


if __name__ == "__main__":
    comm = PMX()
    comm.open('/dev/ttyUSB0', '19200')
    comm.Clear()
    comm2 = PMX()
    comm2.open('/dev/ttyUSB1', '19200')
    comm2.Clear()
    comm3 = PMX()
    comm3.open('/dev/ttyUSB2', '19200')
    comm3.Clear()

    comm.Output("OFF")
    comm2.Output("OFF")
    comm3.Output("OFF")

    comm.close()
    comm2.close()




