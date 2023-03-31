class pyclass:
    def	__init__(self):
        self.x = 0
        self.y = 0

    def set(self, x, y):
        self.x = x
        self.y = y

    def add(self):
        return self.x +	self.y

    def	show(self):
        print( self.x, self.y )

test = pyclass()
test.set(50, -51)
test.show()
print( test.add() )
