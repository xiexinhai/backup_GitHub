# Power supply control interface

**These macros are used for the control of the power supply**

The macros interface_PS_A.py, interface_PS_B.py and interface_PS_C.py
can control a single power supply with command line interface:
```
    print('Welcome to the interface of the power supply!')
    print('Enter 1 to turn on the output.')
    print('Enter 2 to turn off the output.')
    print('Enter 3 to set the voltage value.')
    print('Enter 4 to set the current value.')
    print('Enter 0 to close the interface.')

    while True:
        inNum = input('Please enter the integer number (0~4): \n>')
```

- For close the macro, turn on and turn off the power supply:
```
        if int(inNum) == 0:
            print('Interface macro closed!')
            break

        elif int(inNum) == 1:
            comm.Output("ON")
            print('Power supply output ON!')

        elif int(inNum) == 2:
            print('Power supply output OFF!')
            comm.Output("OFF")
```
- For change the output values of the power supplies:
```
        elif int(inNum) == 3:
            volt = input('Please enter the voltage (0~20V): \n>')
            if float(volt) >= 0 and float(volt)<= 20:
                comm.SetVoltage(float(volt))
                print('Voltage output value set to be U='+str(float(volt))+'V')
            else:
                print('Voltage value error!')

        elif int(inNum) == 4:
            curr = input('Please enter the current (0~5A): \n>')
            if float(curr) >= 0 and float(curr)<= 5:
                comm.SetCurrent(float(curr))
                print('Current output value set to be I='+str(float(curr))+'A')
            else:
                print('Current value error!')
```
The input values for voltage and current is currently limited as $U:[0,20]V$, and $I: [0,5]A$. The limit can be modified.


The macro interface_PS_together.py 
can control three power supplies together.

