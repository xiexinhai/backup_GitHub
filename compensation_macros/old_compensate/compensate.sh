cd ../
python fix_USBID_order.py
/bin/sleep 1
python open_all_PS.py
/bin/sleep 1

while [ true ]; do 
/bin/date
python monitor.py
/bin/sleep 0.2
python control.py
/bin/sleep 1
done

python close_all_PS.py
