# ADC monitor

**These macros are used for monitor of the ADC readout**


The ADC readout will can be saved into a .csv file **log_ADC.csv**, 
by running the shell macro: 

`sh save_ADC.sh`

which contains a simple circulation:
```
while [ true ]; do 
/bin/date
python save_ADC.py
/bin/sleep 1
done
```
It will run the python macro **save_ADC.py** in a time interval as 1 second, which is the saving frequency for ADC readout. 

The python macro **save_ADC.py** will obtain ADC readout and save the values as a list into the .csv file **log_ADC.csv**.
```
def read_aio(aio, channel_id):
    return aio.analog_read_volt(channel_id, aio.DataRate.DR_860SPS)

def read_aio_all(aio):
    ret = []
    for channel in range(32):
        ret.append(read_aio(aio,channel))
    return ret

if __name__ == "__main__":
    aio = aio320ra.AIO_32_0RA_IRC(0x49, 0x3e)
    f = open('log_ADC.csv', "a")
    adc = read_aio_all(aio)
    adc.insert(0,datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"))

    ret = str(adc).strip('[')
    ret = ret.strip(']')
    f.write(ret)
    f.write('\n')
    f.close()
```

The saving order in .csv file is:
`[System time, readout from ADC*[32]]`
