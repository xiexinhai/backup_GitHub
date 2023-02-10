import sys
sys.path.insert(0, '/home/xinhai/MagnetControl/fieldcontrol')

from PMX import PMX
import time
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

    with open("ADCreadout.json","r") as file_obj:
        ADCreadout = json.load(file_obj)
        file_obj.close()
    chx = 7
    chy = 8
    chz = 9


    #magnetic field at the coordinate of censor FM-3500 placed at target center
    Bx = 0.6162593319995725
    By = 0.6547714841263379
    Bz = 0.665111485585729
    #Bzy = -0.0064509222076331715
    #Bzx = -0.006450922207634864

    Ixnow = comm.ReadCurrent()[1]
    Iynow = comm2.ReadCurrent()[1]
    Iznow = comm3.ReadCurrent()[1]

    #Bxnow = Ixnow*Bx*100
    #Bynow = Iynow*By*100
    #Bznow = Iznow*Bz*100

    deltaBx = ADCreadout[chx]*10 #- Bxnow
    Ix = Ixnow + deltaBx/Bx/100

    deltaBy = ADCreadout[chy]*10 #- Bynow
    Iy = Iynow + deltaBy/By/100

    deltaBz = ADCreadout[chz]*10 #- Bznow
    Iz = Iznow + deltaBz/Bz/100


    #Iz = ADCreadout[chz]/Bz/10

    #Iy = (ADCreadout[chy]*10 - Iz*Bzy*100)/By/100

    #Ix = (ADCreadout[chx]*10 - Iz*Bzx*100)/Bx/100


#    print("signal =",ADCreadout[7],"Gs")
#    print("Ix =",Ix,"A")
#    print("Iy =",Iy,"A")
#    print("Iz =",Iz,"A")

#    comm.SetVoltage(Ix)
#    comm2.SetVoltage(Iy)
#    comm3.SetVoltage(Iz)

    comm.SetCurrent(Ix)
    comm2.SetCurrent(Iy)
    comm3.SetCurrent(Iz)

    comm.close()
    comm2.close()
    comm3.close()



