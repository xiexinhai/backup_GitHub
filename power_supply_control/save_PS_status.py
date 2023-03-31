import sys
import os
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

    f = open('./current_log/log_PS.csv', "a")
    f.write('time,V1,I1,V2,I2,V3,I3')
    f.write('\n')
    f.close()

    count = 0
    if_title = False

    try:
        while True:
            f = open('./current_log/log_PS.csv', "a")
            if if_title == True:
                f.write('time,V1,I1,V2,I2,V3,I3')
                f.write('\n')
                if_title = False

            stat = []
        
            stat.append(comm.ReadVoltage()[1])
            stat.append(comm.ReadCurrent()[1])
            stat.append(comm2.ReadVoltage()[1])
            stat.append(comm2.ReadCurrent()[1])
            stat.append(comm3.ReadVoltage()[1])
            stat.append(comm3.ReadCurrent()[1])
        
            #stat.insert(0,datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))
            stat.insert(0,datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            print(stat[0])
        
            ret = str(stat).strip('[')
            ret = ret.strip(']')
            f.write(ret)
            f.write('\n')
            f.close()
            count += 1

            if count >= 1800:
                os.rename("./current_log/log_PS.csv","saved_log/log_PS_"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+".csv")
                count = 0
                if_title = True

            time.sleep(2)

    except KeyboardInterrupt:
        print("Stopped by KeyboardInterrupt!")
        f.close()
        os.rename("./current_log/log_PS.csv","saved_log/log_PS_"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+".csv")
    
    except Exception as e:
        print(e)
        f.close()
        os.rename("./current_log/log_PS.csv","saved_log/log_PS_"+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+".csv")

    comm.close()
    comm2.close()
    comm3.close()



