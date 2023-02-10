# Compensation system

**The macros are used to compensate the magnetic field**

Current setup is: Put a magnetic sensor in the center of the chamber, and the comensate system will generate magnetic field to compensate the three probes of the sensor to zero magnetic field.
**The setup need to be modified under different situations**.

One can run the macro by:
`python control.py`

**Pre-definition of the compensation macro:**
- Define the time interval and ADC channel numbers for three-axis probes and temperature sensor:
```
deltaT = 5.0

Bchx = 1 # LR coil
Bchy = 0 # UD coil
Bchz = 2 # FB coil

T1ch = 14
T2ch = 15
```
- Define function to obtain ADC readout:
```
def read_aio(aio, channel_id):
    return aio.analog_read_volt(channel_id, aio.DataRate.DR_860SPS)

def read_aio_all(aio):
    ret = []
    for channel in range(32):
        ret.append(read_aio(aio,channel))
    return ret
```
- Define funtion for safety check:
```
def check_safety(ADCreadout,Ix,Iy,Iz):
    ret = True

    #Current check
    if Ix > 5. or Iy > 5. or Iz > 5.:
        ret = False

    #Temperature check
    T1 = abs(ADCreadout[T1ch] - ADCreadout[T1ch+16])
    T2 = abs(ADCreadout[T2ch] - ADCreadout[T2ch+16])
    if T1 > 30. or T2 > 30.:
        ret = False

    #Probe signal check (FM-3500)
    Bx = 10.*abs(ADCreadout[Bchx] - ADCreadout[Bchx+16])
    By = 10.*abs(ADCreadout[Bchy] - ADCreadout[Bchy+16])
    Bz = 10.*abs(ADCreadout[Bchz] - ADCreadout[Bchz+16])
    if Bx > 100 or By > 100 or Bz > 100:
        ret = False

    return ret
```
The safety check will be False under such situations: Needed output current larger than $5A$, temperature is higher than 30 degree, probes detected abnormal large magnetic field, larger than 100 $\mu T$. Other cases, the safety check will be True.

- Define function to get $\Delta B$:
```
def get_deltaB(ADCreadout):
    #need to be modified when upgrading probes!!!

    deltaBx = 10*(ADCreadout[Bchx] - ADCreadout[Bchx+16])
    deltaBy = -1*10*(ADCreadout[Bchy] - ADCreadout[Bchy+16])
    deltaBz = 10*(ADCreadout[Bchz] - ADCreadout[Bchz+16])

    return (deltaBx, deltaBy, deltaBz)
```
Under current setup, the $\Delta B$ values are directly the signals from probes.

- Deine function to calculate the needed compensation current based on $\Delta B$ for three coils:
```
def cal_deltaI(deltaBx,deltaBy,deltaBz):
    #need to be modified when selecting different compensation coordinate!!!

    #magnetic field at the coordinate of censor FM-3500 placed at target center
    #x-probe (0,0.014,-0.0315)
    #y-probe (0.014,0,-0.0035)
    #z-probe (-0.014,-0.014,-0.0175)

    #initial: Bx0 ~ 27, By0 ~ -25, Bz0 ~ -5

    Bx  = -1 * 0.62432752
    Bxz = -1 * -0.003345779

    By  = 0.657061546
    Byz = -0.003588381

    Bz  = -1 * -0.65464895

    deltaIx = -deltaBx/Bx/100

    deltaIy = -deltaBy/By/100

    deltaIz = -deltaBz/Bz/100 - deltaIx*Bxz/Bz - deltaIy*Byz/Bz

    return (deltaIx,deltaIy,deltaIz)
```

**Main structure of the compensation macro:**
- Open Serial communication of the power supplies according to the order UD, LR and FB, set the over U&I protection:
```
if __name__ == "__main__":
    #When rebot the server, the USB ID of the power supply may change!!!!!!
    with open("../USBID_order.json","r") as file_obj:
        USBID_order = json.load(file_obj)
        file_obj.close()

    #UD coil
    comm = PMX()
    comm.open(USBID_order[0], '19200')
    comm.Clear()
    comm.SetOverCurrentProtection(5.0)
    comm.SetOverVoltageProtection(20)
    
    #LR coil
    comm2 = PMX()
    comm2.open(USBID_order[1], '19200')
    comm2.Clear()
    comm2.SetOverCurrentProtection(5.0)
    comm2.SetOverVoltageProtection(20)

    #FB coil
    comm3 = PMX()
    comm3.open(USBID_order[2], '19200')
    comm3.Clear()
    comm3.SetOverCurrentProtection(5.0)
    comm3.SetOverVoltageProtection(20)
```

- The main circulation structure:
```
    try:
        comm.SetCurrent(0.)
        comm.Output("ON")
        comm2.SetCurrent(0.)
        comm2.Output("ON")
        comm3.SetCurrent(0.)
        comm3.Output("ON")
        while True:
            #ADC readout
            aio = aio320ra.AIO_32_0RA_IRC(0x49, 0x3e)
            ADCreadout = read_aio_all(aio)
        
            #current values for power supply
            Ixnow = comm.ReadCurrent()[1]
            Iynow = comm2.ReadCurrent()[1]
            Iznow = comm3.ReadCurrent()[1]
        
            deltaBx, deltaBy, deltaBz = get_deltaB(ADCreadout)

            deltaIx, deltaIy, deltaIz = cal_deltaI(deltaBx,deltaBy,deltaBz)

            Ix = Ixnow + deltaIx
            Iy = Iynow + deltaIy
            Iz = Iznow + deltaIz

            if check_safety(ADCreadout,Ix,Iy,Iz) == False:
                break

            comm.SetCurrent(Iy)
            comm2.SetCurrent(Ix)
            comm3.SetCurrent(Iz)

            print("--------------------Ctrl+C to stop------------------------")
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))
            print("Ix =",Ix,"A")
            print("Iy =",Iy,"A")
            print("Iz =",Iz,"A")

            #compensation frequency/time-interval!!!
            time.sleep(deltaT)
```
The process is: 
1. Turn on three power supplies
2. Enter the circulation
3. Obtain ADC readout and current output values of power supplies
4. Obtain $\Delta B$ from ADC readout
5. Calculate the needed output current for three power supplies
6. Safety check
7. (A) Safety check is True, have a system time sleep (e.g. $1s$, can be modified), then go back to process 3.
7. (B) Safety check is False, circulation will quit.

- When have KeyboardInterrupt or macro crash, turn off all power supplies:
```
    except KeyboardInterrupt:
        print("Stopped by KeyboardInterrupt!")
        comm.Output("OFF")
        comm2.Output("OFF")
        comm3.Output("OFF")

        comm.SetCurrent(0.)
        comm2.SetCurrent(0.)
        comm3.SetCurrent(0.)

    except Exception as e:
        print(e)
        comm.Output("OFF")
        comm2.Output("OFF")
        comm3.Output("OFF")

        comm.SetCurrent(0.)
        comm2.SetCurrent(0.)
        comm3.SetCurrent(0.)


    print("Stopped by safety check!")
    comm.Output("OFF")
    comm2.Output("OFF")
    comm3.Output("OFF")

    comm.SetCurrent(0.)
    comm2.SetCurrent(0.)
    comm3.SetCurrent(0.)

    comm.close()
    comm2.close()
    comm3.close()
```





