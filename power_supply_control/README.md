# Power supply control

**print_PS_status.py:** Print out once the current power supplies information

**save_PS_status.py:** Saving the power supplies information into .csv file

**command_control.py:**
Power supply control using command line

Select which coil:
```
python command_control.py -c UD 
python command_control.py -c LR
python command_control.py -c FB
```

Switch coil to turn ON or turn OFF: (e.g. turn ON or turn OFF the UD coil)
```
python command_control.py -c UD -s ON
python command_control.py -c UD -s OFF
```

Change the current output: (e.g. change the UD power supply output current to 0.2 A)
```
python command_control.py -c UD -m current -i 0.2
```

Change the current output: (e.g. change the UD power supply output voltage to 0.2 V)
```
python command_control.py -c UD -m voltage -v 0.2
```
