# File:"C:\Users\DELL\Desktop\PhD\energy function\R1\gaus_siedel_time1.py", generated on WED, SEP 30 2015  23:28, release 33.05.02
psspy.read(0,r"""C:\Users\DELL\Desktop\PhD\energy function\NE\ieee39bus.raw""")
psspy.fnsl([0,0,0,1,0,0,99,0])
psspy.branch_chng(1,2,r"""1""",[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])



import time
time_start = time.clock()

psspy.fnsl([0,0,0,1,0,0,99,0])

time_elapsed = (time.clock() - time_start)
print time_elapsed



#psspy.mslv([0,0,0,1,0,0])
#psspy.solution_parameters_4([_i,_i,_i,_i,10],[_f,_f,_f, 0.0001,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])

#psspy.branch_chng(1,2,r"""1""",[1,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
#psspy.solv([0,0,0,1,0,0])
