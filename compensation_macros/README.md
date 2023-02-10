# Compensate macros
**These macros are used for the compensation system in the g-2/EDM experiment**

The relevent macros can control the power supplies of three coils:
- Up-down (UD) coil
- Left-right (LR) coil
- Forward-backward (FB) coil

Before using the macros of the compensation system,
one need to firstly run the macro: 

`sh fix_USBID_order.py`

which will read the version ID for three power supplies and fix the order into file, **USBID_order.json**, which contains a list, for example:

`["/dev/ttyUSB0", "/dev/ttyUSB1", "/dev/ttyUSB2"]`

For other macros, when it need to use Serial interface to control the power supply, this order will be used, as the order as **UD, LR, FB**, respectively.

The macro **fix_USBID_order.py** is used to save order of three power supplies:

- Define the order of power supply by the veision ID, and construct empty list:
(Here, the string "USB0" need to be finally changed the real Version ID of the power supply for UD coil, and similar for others)
```
    PMXorder = ["USB0","USB1","USB2"]
    fixed_order = []
```

- Use Serial interface to open the communication to the three power supplies, and check each Version ID:
```
    comm = PMX()
    comm.open('/dev/ttyUSB0', '19200')
    comm.Clear()
    ID1 = comm.CheckVersion()[1]

    comm2 = PMX()
    comm2.open('/dev/ttyUSB1', '19200')
    comm2.Clear()
    ID2 = comm2.CheckVersion()[1]

    comm3 = PMX()
    comm3.open('/dev/ttyUSB2', '19200')
    comm3.Clear()
    ID3 = comm3.CheckVersion()[1]
```

- Compare the Version ID with the order we need and fix the order into the previous empty list:

```
    if(ID1 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[0]):
        fixed_order.append("/dev/ttyUSB2")

    if(ID1 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[1]):
        fixed_order.append("/dev/ttyUSB2")

    if(ID1 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB0")
    if(ID2 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB1")
    if(ID3 == PMXorder[2]):
        fixed_order.append("/dev/ttyUSB2")
```

- Save the fixed order as a list into the .json file and close the Serial interface:
```
    with open("USBID_order.json", 'w') as file_obj:
        json.dump(fixed_order, file_obj)
        file_obj.close()

    comm.close()
    comm2.close()
    comm3.close()
```












