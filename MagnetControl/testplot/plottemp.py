# Library
import numpy as np
import matplotlib.pyplot as plt
import collections

# params
dataLength = 1000  # １つのデータの配列の点数
frame = 1000  # プロットするフレーム数
sleepTime = 1.0  # １フレーム表示する時間[s]

history = collections.deque(maxlen=dataLength)

# plot
if __name__== "__main__":
    try:
        for i in range(frame): # フレーム回数分グラフを更新
            f = open('/sys/class/thermal/thermal_zone0/temp','r')
            data = float( f.read() )# プロットするデータを作成
            print(data)
            f.close()
            
            history.append(data)
            x = list(range(i-len(history), i))
            plt.plot(x,history, label='channel 0') # データをプロット
            plt.legend()
            plt.xlabel("second")
            plt.ylabel("Temperature [Deg.C]")
            plt.draw() # グラフを画面に表示開始
            plt.pause(sleepTime) # SleepTime時間だけ表示を継続
            plt.cla() # プロットした点を消してグラフを初期化
    except KeyboardInterrupt:
        plt.close()
plt.close()