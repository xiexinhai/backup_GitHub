import serialcomm.py
import time

# Open and initialize the serial communication to power supplies
pmx_x = SerialComm()
pmx_x.open('/dev/ttyUSB0', '19200')

pmx_x.send('IDN?')
result, data = pmx_x.recv(10)
print(result)
print(data)

#Set voltage
Vx  = 0.0
Vy  = 0.0
Vz  = 0.0
pmx_x.SetVoltage(Vx)
pmx_y.SetVoltage(Vy)
pmx_z.SetVoltage(Vz)

#Set target field
B0x = 0.
B0y = 0.
B0z = 0.

#Set default voltage to each coil

while True:
    #Measure B-field
    atoB = 1.0
    
    adc = ADC.read()

    Bx  = atoB * ( adc[ix] - ped[ix] )
    By  = atoB * ( adc[iy] - ped[iy] )
    Bz  = atoB * ( adc[iz] - ped[iz] )

    dBx = Bx - B0x
    dBy = By - B0y
    dBz = Bz - B0z

    dVx = dBx / atoB
    dVy = dBy / atoB
    dVz = dBz / atoB

    #Read current voltage
    V0x = pmx_x.Voltage()
    V0y = pmx_y.Voltage()
    V0z = pmx_z.Voltage()

    Vx  = V0x + dVx
    Vy  = V0y + dVy
    Vz  = V0z + dVz
    
    #Apply corrections to voltage
    pmx_x.SetVoltage(Vx)
    pmx_y.SetVoltage(Vy)
    pmx_z.SetVoltage(Vz)

    time.sleep(5)
    
    if IsStop:
        break

