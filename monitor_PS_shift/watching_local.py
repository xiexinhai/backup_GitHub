import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time


if __name__ == "__main__":
    ct = None
    count_ct = 0

    while True:
        try:
            df = pd.read_csv("log/log_PS.csv")
            t = df["time"]
            t = np.array(t)
        
            arr = []
            arr.append(np.array(df['V1']))
            arr.append(np.array(df['I1']))
            arr.append(np.array(df['V2']))
            arr.append(np.array(df['I2']))
            arr.append(np.array(df['V3']))
            arr.append(np.array(df['I3']))

            for i in range(len(arr)):
                np.delete(arr[i],-1)
    
            print('--------------------------------')
            print(t[-1])
            print('UD power supply:')
            print('Voltage output = '+str(arr[0][-2])+' V')
            print('Current output = '+str(arr[1][-2])+' A')
            print('')
            print('LR power supply:')
            print('Voltage output = '+str(arr[2][-2])+' V')
            print('Current output = '+str(arr[3][-2])+' A')
            print('')
            print('FB power supply:')
            print('Voltage output = '+str(arr[4][-2])+' V')
            print('Current output = '+str(arr[5][-2])+' A')
            print('--------------------------------')

            n_lim = 120
            I1 = []
            I2 = []
            I3 = []
            t_ls = []
            if len(t) > n_lim+5:
                for i in range(n_lim):
                    temp = str(t[-n_lim+i-2])
                    temp = temp.replace("'","")
                    t_ls.append(str(temp))
                    I1.append(arr[1][-n_lim+i-2])
                    I2.append(arr[3][-n_lim+i-2])
                    I3.append(arr[5][-n_lim+i-2])

            else:
                for i in range(len(t)-2):
                    temp = str(t[i])
                    temp = temp.replace("'","")
                    t_ls.append(str(temp))
                    I1.append(arr[1][i])
                    I2.append(arr[3][i])
                    I3.append(arr[5][i])

            tmp_ct = t_ls[-1]
            if tmp_ct == ct:
                if count_ct < 10:
                    count_ct += 1
            else:
                count_ct = 0

            ct = tmp_ct
            if count_ct >= 10:
                print('========================================')
                print('========================================')
                print('Log file WARNING!!')
                print('Log file not updated for more than 30s!!!')
                print('Check the log saving macro in the Raspberry Pi!!!')
                print('========================================')
                print('========================================')

            maxI1 = max(I1)
            idMaxI1 = [i for i, j in enumerate(I1) if j == maxI1][0]
            maxI2 = max(I2)
            idMaxI2 = [i for i, j in enumerate(I2) if j == maxI2][0]
            maxI3 = max(I3)
            idMaxI3 = [i for i, j in enumerate(I3) if j == maxI3][0]

            minI1 = min(I1)
            idMinI1 = [i for i, j in enumerate(I1) if j == minI1][0]
            minI2 = min(I2)
            idMinI2 = [i for i, j in enumerate(I2) if j == minI2][0]
            minI3 = min(I3)
            idMinI3 = [i for i, j in enumerate(I3) if j == minI3][0]

            th = 0.1
            if abs(maxI1-minI1) >= th:
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('Up-down coil power supply WARNING!')
                print('Maximum output current: '+str(maxI1)+' A, time: '+str(t_ls[idMaxI1]))
                print('Minimum output current: '+str(minI1)+' A, time: '+str(t_ls[idMinI1]))
                print('Variation of current output larger than '+str(th)+' A!!!!')
                print('Time period: from '+str(t_ls[0])+' to '+str(t_ls[-1]))
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')

            if abs(maxI2-minI2) >= th:
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('Up-down coil power supply WARNING!')
                print('Maximum output current: '+str(maxI2)+' A, time: '+str(t_ls[idMaxI2]))
                print('Minimum output current: '+str(minI2)+' A, time: '+str(t_ls[idMinI2]))
                print('Variation of current output larger than '+str(th)+' A!!!!')
                print('Time period: from '+str(t_ls[0])+' to '+str(t_ls[-1]))
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')

            if abs(maxI3-minI3) >= th:
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('Up-down coil power supply WARNING!')
                print('Maximum output current: '+str(maxI3)+' A, time: '+str(t_ls[idMaxI3]))
                print('Minimum output current: '+str(minI3)+' A, time: '+str(t_ls[idMinI3]))
                print('Variation of current output larger than '+str(th)+' A!!!!')
                print('Time period: from '+str(t_ls[0])+' to '+str(t_ls[-1]))
                print('**************************************************************************************')
                print('**************************************************************************************')
                print('**************************************************************************************')

            time.sleep(3)

        except KeyboardInterrupt:
            print("Stopped by KeyboardInterrupt!")
            break
    
        except Exception as e:
            #print(e)
            continue

