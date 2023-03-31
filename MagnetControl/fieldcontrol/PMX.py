import time
import struct
from serial.tools import list_ports
from SerialComm import SerialComm    # SerialComm

class PMX(SerialComm):
    
    def __init__(self):
        super().__init__()
        
    def Clear(self):
        command = '*CLS'
        super().send(command)
        return
    
    def CheckVersion(self):
        command = '*IDN?'
        super().send(command)
        result,data = super().recv(1)
        ddata = data.strip().decode('UTF-8')
        return (result,ddata)

    def SetVoltage(self, V):
        command = 'VOLT ' + str(V)
        super().send(command)
        return

    def SetCurrent(self, I):
        command = 'CURR ' + str(I)
        super().send(command)
        return

    def SetOverCurrentProtection(self, I):
        command = 'CURR:PROT ' + str(I)
        super().send(command)
        return

    def SetOverVoltageProtection(self, V):
        command = 'CURR:PROT ' + str(V)
        super().send(command)
        return
    
    def ReadVoltage(self):
        command = 'MEAS:VOLT?'
        super().send(command)
        result,data = super().recv(1)
        ddata = float( data.strip().decode('UTF-8') )
#        print(ddata, type(ddata))
        return (result,ddata)
    
    def ReadCurrent(self):
        command = 'MEAS:CURR?'
        super().send(command)
        result,data = super().recv(1)
        ddata = float( data.strip().decode('UTF-8') )
#        print(ddata, type(ddata))
        return (result,ddata)
    
    def Output(self, IsON):
        command = 'OUTP ' + IsON
        super().send(command)
        return

    def CheckError(self):
        command = 'SYST:ERR?'
        super().send(command)
        result,data = super().recv(1)
        ddata = data.strip().decode('UTF-8')
        return (result,ddata)



if __name__ == "__main__":

    comm = PMX()
    comm.open('/dev/ttyUSB0', '19200')
    comm.Clear()

    result, data = comm.CheckVersion()
    print("ID?")
    print(result, data)

    comm.SetVoltage(3.05)
    time.sleep(1)
    result, V = comm.ReadVoltage()
    result, I = comm.ReadCurrent()
    print("Vmon=", V, " Imon=", I )
    
    time.sleep(1)
    comm.Output("ON")    #    コマンド送信
    time.sleep(1)
    result, V = comm.ReadVoltage()
    result, I = comm.ReadCurrent()
    print("Vmon=", V, " Imon=", I )
    
    time.sleep(3)
    comm.Output("OFF")    #    コマンド送信

    time.sleep(3)
    result, V = comm.ReadVoltage()
    result, I = comm.ReadCurrent()
    print("Vmon=", V, " Imon=", I )

    time.sleep(3)
    result, V = comm.ReadVoltage()
    result, I = comm.ReadCurrent()
    print("Vmon=", V, " Imon=", I )

    result, error = comm.CheckError()
    print("Status: ", error )
    
    # シリアルを閉じる
    comm.close();
