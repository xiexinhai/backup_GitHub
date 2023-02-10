# Calculation of magnetic field

**These macros are used for calculation of magnetic field of the coils, UD, LR and FB**

The macros are divied into three .py parts, stand for three coils:
- B_up_down.py
- B_left_right.py
- B_forward_back.py

Give a position, give a current, and the macro can calculate the generated magnetic field at the given position. 

An example:
```
from B_left_right import B_left_right
from B_forward_back import B_forward_back
from B_up_down import B_up_down
import math

if __name__ == "__main__":
    a = B_left_right()
    a.set_I(1)
    a.set_xyz(0.0,0.0,0.0)
    print(a.magnetic3D())

    b = B_up_down()
    b.set_I(1)
    b.set_xyz(0.0,0.0,0.0)
    print(b.magnetic3D())

    c = B_forward_back()
    c.set_I(1)
    c.set_xyz(0.0,0.0,0.0)
    print(c.magnetic3D())
```
This will calculate the the generate magnetic field at the center $(0,0,0)$ with the current $I=1.0A$.