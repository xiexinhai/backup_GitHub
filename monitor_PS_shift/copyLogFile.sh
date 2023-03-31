#!/usr/bin/expect -f
set password shift
spawn scp -r shift@10.105.50.62:/home/shift/power_supply_control/current_log/log_PS.csv ./log
expect "*password*"
send "$password\r"
expect eof
