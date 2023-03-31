import sys
sys.path.insert(0, '/home/shift/MagnetControl/fieldcontrol')
import json

from PMX import PMX


if __name__ == "__main__":
    #Version ID for the power supplies
    PMXorder = ["KIKUSUI,PMX18-5A,CQ003823,IFC01.53.0012 IOC01.10.0070","KIKUSUI,PMX18-5A,CQ003824,IFC01.53.0012 IOC01.10.0070","KIKUSUI,PMX18-5A,CQ003961,IFC01.53.0012 IOC01.10.0070"]
    fixed_order = []

    comm = PMX()
    comm.open('/dev/ttyUSB0', '19200')
    comm.Clear()
    ID1 = comm.CheckVersion()[1]

    comm2 = PMX()
    comm2.open('/dev/ttyUSB1', '19200')
    comm2.Clear()
    ID2 = comm2.CheckVersion()[1]

    comm3 = PMX()
    comm3.open('/dev/ttyUSB2', '19200')
    comm3.Clear()
    ID3 = comm3.CheckVersion()[1]

    if(ID1 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB2")

    if(ID1 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB2")

    if(ID1 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB2")

    with open("USBID_order.json", 'w') as file_obj:
        json.dump(fixed_order, file_obj)
        file_obj.close()

    comm.close()
    comm2.close()
    comm3.close()




