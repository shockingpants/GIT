import cellsim as cs
import math
space=cs.cellspace()
space.add_cell((200,100))
space.add_cell((200,250),angle=math.pi/2)
space.run()
#space.stop()
