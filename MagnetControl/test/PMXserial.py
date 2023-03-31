import serial    #    PySerial
from serial.tools import list_ports
import time

dev = "/dev/ttyUSB0"    #    デバイス名
rate = 19200    #    レート (bps)
ser = serial.Serial(dev, rate, timeout=1)

#string = "OUTPut 0"
#string = "OUTPut?"
#string = "MEASure:CURRent?"
string = "MEASure:VOLTage?"
#string = "*IDN?"
string = string + "\n"    #    ターミネーターを付ける
print (string)
#ser.write(str.encode("*IDN?\n"))    #    コマンド送信
ser.write(str.encode(string))    #    コマンド送信

res = ser.readline()    #    コマンド受信
res_disp = res.strip().decode() 
print (res)
print (res_disp)

ser.close()
