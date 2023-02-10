# Power supply monitor

**These macros are used for monitor of the power supplies**

The output information of three power supplies can be saved into a .csv file **log_PS_status.csv**, 
by running the shell macro: 

`sh save_PS_status.sh`

which contains a simple circulation:
```
while [ true ]; do 
/bin/date
python save_PS_status.py
/bin/sleep 1
done
```
It will run the python macro **save_PS_status.py** in a time interval as 1 second, which is the saving frequency for ADC readout. (Time interval can be changed.)

The python macro **save_PS_status.py** will obtain information of the power supplies and save the values as a list into the .csv file **log_PS_status.csv**:

- Read the previous fixed order for three power supplies
```
    with open("../USBID_order.json","r") as file_obj:
        USBID_order = json.load(file_obj)
        file_obj.close()
```
- Open the Serial communication for power supplies, open .csv file, and creat an empty list
```
    comm = PMX()
    comm.open(USBID_order[0], '19200')
    comm.Clear()
    
    comm2 = PMX()
    comm2.open(USBID_order[1], '19200')
    comm2.Clear()

    comm3 = PMX()
    comm3.open(USBID_order[2], '19200')
    comm3.Clear()

    f = open('log_PS_status.csv', "a")
    stat = []
```
- Save the output information and close Serial communication
```
    stat.append(comm.ReadVoltage()[1])
    stat.append(comm.ReadCurrent()[1])
    stat.append(comm2.ReadVoltage()[1])
    stat.append(comm2.ReadCurrent()[1])
    stat.append(comm3.ReadVoltage()[1])
    stat.append(comm3.ReadCurrent()[1])

    stat.insert(0,datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))

    ret = str(stat).strip('[')
    ret = ret.strip(']')

    comm.close()
    comm2.close()
    comm3.close()
```

The saving order in .csv file is:
`[System time, U(PS_UD), I(PS_UD), U(PS_LR), I(PS_LR), U(PS_FB), I(PS_FB)]`

If one just need to print out the information of power supplies once, one can run the macro:

`python print_PS_status.py`



