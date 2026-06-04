# File:"C:\Users\DELL\Desktop\PhD\energy function\R1\gaus_siedel_time.py", generated on WED, SEP 30 2015  22:59, release 33.05.02

psspy.read(0,r"""C:\Users\DELL\Desktop\PhD\energy function\NE\ieee39bus.raw""")

import time
time_start = time.clock()

psspy.fnsl([0,0,0,1,0,0,99,0])

time_elapsed = (time.clock() - time_start)
print time_elapsed

#psspy.branch_chng(1,2,r"""2""",[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
#psspy.solv([0,0,0,0,0,0])




#psspy.branch_chng(1,2,r"""2""",[1,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
#psspy.mslv([0,0,0,0,0,0])
