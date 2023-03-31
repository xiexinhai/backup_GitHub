import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import sys

if __name__ == "__main__":

    while True:
        try:
            #os.system('cp log/log_PS.csv log/temp.csv')
            #df = pd.read_csv("log/temp.csv")
            df = pd.read_csv("log/log_PS.csv")
            t = df["time"]
            t = np.array(t)
            t_ls = []
        
            arr = []
            arr.append(np.array(df['V1']))
            arr.append(np.array(df['I1']))
            arr.append(np.array(df['V2']))
            arr.append(np.array(df['I2']))
            arr.append(np.array(df['V3']))
            arr.append(np.array(df['I3']))
        
            n_lim = 300
            V1=[]
            I1=[]
            V2=[]
            I2=[]
            V3=[]
            I3=[]
            if len(t) > n_lim+5:
                for i in range(n_lim):
                    temp = str(t[-n_lim+i-2])
                    temp = temp.replace("'","")
                    t_ls.append(str(temp))
                    V1.append(arr[0][-n_lim+i-2])
                    I1.append(arr[1][-n_lim+i-2])
                    V2.append(arr[2][-n_lim+i-2])
                    I2.append(arr[3][-n_lim+i-2])
                    V3.append(arr[4][-n_lim+i-2])
                    I3.append(arr[5][-n_lim+i-2])
                    
            else:
                for i in range(len(t)-2):
                    temp = str(t[i])
                    temp = temp.replace("'","")
                    t_ls.append(str(temp))
                    V1.append(arr[0][i])
                    I1.append(arr[1][i])
                    V2.append(arr[2][i])
                    I2.append(arr[3][i])
                    V3.append(arr[4][i])
                    I3.append(arr[5][i])
    
            plt.close()
            fig = plt.figure(figsize=(10, 7.5))
            fig.subplots_adjust(wspace=0.2,hspace=0.5)
            subfig = []
    
            #for i in range(6):
            #    subfig.append(fig.add_subplot(6,1,i+1))
            #    subfig[i].set_xticks([])
    
            #subfig[0].set_ylabel('U_UD(V)')
            #subfig[1].set_ylabel('I_UD(A)')
            #subfig[2].set_ylabel('U_LR(V)')
            #subfig[3].set_ylabel('I_LR(A)')
            #subfig[4].set_ylabel('U_FB(V)')
            #subfig[5].set_ylabel('I_FB(A)')
    
            subfig.append(fig.add_subplot(2,1,1))
            subfig[0].set_xticks([])
            subfig.append(fig.add_subplot(2,1,2))
            subfig[1].set_xticks([])
        
            subfig[0].set_ylabel('U(V)')
            subfig[1].set_ylabel('I(A)')
    
            tick = 1
            if len(t_ls)<10:
                tick = 1
            else:
                tick = int(len(t_ls)/10)
    
            xlabel = [i for i in range(0,len(t_ls),tick)]
            if (len(t_ls)-1 in xlabel) == False:
                xlabel.append(len(t_ls)-1)
        
            #subfig[0].plot(V1, markersize=4,marker='.',label='V1',c='r')
            #subfig[1].plot(I1, markersize=4,marker='.',label='I1',c='r')
            #subfig[2].plot(V2, markersize=4,marker='.',label='V2',c='r')
            #subfig[3].plot(I2, markersize=4,marker='.',label='I2',c='r')
            #subfig[4].plot(V3, markersize=4,marker='.',label='V3',c='r')
            #subfig[5].plot(I3, markersize=4,marker='.',label='I3',c='r')
    
            subfig[0].plot(V1, markersize=4,marker='.',label='UD',c='r')
            subfig[0].plot(V2, markersize=4,marker='*',label='LR',c='b')
            subfig[0].plot(V3, markersize=4,marker='x',label='FB',c='k')
            subfig[0].legend()
    
            subfig[1].plot(I1, markersize=4,marker='.',label='UD',c='r')
            subfig[1].plot(I2, markersize=4,marker='*',label='LR',c='b')
            subfig[1].plot(I3, markersize=4,marker='x',label='FB',c='k')
            subfig[1].legend()
    
            plt.xticks(xlabel,[t_ls[i] for i in xlabel], rotation=90)
            plt.tight_layout()
    
            #thread1 = Thread(target=close, args=(3,))
            plt.pause(10)
            #plt.close()

        except KeyboardInterrupt:
            print("Stopped by KeyboardInterrupt!")
            plt.close()
            break
    
        except Exception as e:
            #print(e)
            continue

