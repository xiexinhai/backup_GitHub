# Run fist to fix ID for three power supplies
**Run fist to fix ID for three power supplies**

The relevent macros can control the power supplies of three coils:
- Up-down (UD) coil
- Left-right (LR) coil
- Forward-backward (FB) coil

Before using the macros of the compensation system,
one need to firstly run the macro: 

`python fix_USBID_order.py`

which will read the version ID for three power supplies and fix the order into file, **USBID_order.json**, which contains a list, for example:

`["/dev/ttyUSB0", "/dev/ttyUSB1", "/dev/ttyUSB2"]`

For other macros, when we need to use Serial interface to control the power supply, this order will be used, as the order as **UD, LR, FB**, respectively.


