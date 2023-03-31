# Power supply monitoring
**These macros are used to monitor the power supplies for three coils**

These macros can run at local computer instead of Raspberry Pi.

**run_copyLog.sh:** Copy log file from Raspberry Pi to the local (time interval 1s, also depend on network speed).
```
sh run_copyLog.sh
```

**plot_local.py:** Drawing plots of the information for three power supplies. Updating time interval 5s, also depend on network speed.
```
python plot_local.py
```

**watching_local.py:** Spectators, print out the recent log, and can give warnings under certain situation.

