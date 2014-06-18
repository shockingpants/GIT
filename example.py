import cellsim as cs
import math
space=cs.cellspace()
space.add_cell((200,300))
#space.add_cell((200,240),angle=math.pi/2)
space.run()
#space.stop()
