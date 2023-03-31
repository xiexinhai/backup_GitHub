import PMX

pmx1 = PMX()
pmx1.send("IDN?")
print( pmx1.read() )
