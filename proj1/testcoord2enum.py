#!/usr/bin/env python3

from cyl2d import coord2enum, enum2coord

for i in range(100):
    for j in range(100):

        enum = coord2enum(i,j,100)
        coord = enum2coord(enum, 100)

        try:

            assert coord == (i,j)
        except:
            print("mismatch")
            print('should be',i,j)
            print('enum is',enum)
            print('actually got',coord)
            quit()
