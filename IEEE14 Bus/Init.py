# File:"C:\Users\Yang\Documents\PTI\PSSE34\test-IEEE\IEEE9busSystem - fix bug\Init.py", generated on WED, MAY 16 2018  15:53, PSS(R)E release 34.03.02

import psse3603  # type: ignore
import psspy     # type: ignore
import os
import sys
psspy.psseinit(50)
raw_file = r"""C:\Users\Acer\Desktop\IEEE14 Bus\IEEE14bus_v32.raw"""
dyr_file = r"""C:\Users\Acer\Desktop\IEEE14 Bus\ieee14.dyr"""
out_file= r"""C:\Users\Acer\Desktop\IEEE14 Bus\ieee14.out"""
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.fact()
psspy.tysl(0)
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,12,0])
psspy.chsb(0,1,[-1,-1,-1,1,13,0])
psspy.chsb(0,1,[-1,-1,-1,1,16,0])

#channel Setup
psspy.delete_all_plot_channels()

