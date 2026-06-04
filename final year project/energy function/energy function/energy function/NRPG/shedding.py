# File:"C:\Users\DELL\Desktop\PhD\energy function\NRPG\shedding.py", generated on WED, OCT 07 2015  19:50, release 32.00.03
P=360
Xd=0.2297
H=3.7
shed=0.80
remain=1-shed
Premain=P*remain
Pshed=P*shed
Xdremain=Xd*(1+shed/remain)
Xdshed=Xd*(1+remain/shed)
Hremain=H*remain
Hshed=H*shed
Pmax=588.
Qmin=-36.88
tf=1.39
delay=tf+0.0+2*0.0

psspy.case(r"""C:\Users\DELL\Desktop\PhD\energy function\NRPG\NRexciter_pristine_cls.sav""")
psspy.machine_data_2(22,r"""2""",[_i,_i,_i,_i,_i,_i],[ Pshed,_f,_f,_f,_f,_f, Pmax,_f, Xdshed,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_data_2(22,r"""1""",[_i,_i,_i,_i,_i,_i],[ Premain,Qmin,_f,_f,_f,_f,_f,_f, Xdremain,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.dyre_new([1,1,1,1],r"""C:\Users\DELL\Desktop\PhD\energy function\NRPG\NRexciter_pristine_cls.dyr""","","","")
psspy.add_plant_model(22,r"""2""",1,r"""GENCLS""",0,"",0,[],[],2,[0.0,0.0])
psspy.change_plmod_con(22,r"""2""",r"""GENCLS""",1, Hremain)
psspy.change_plmod_con(22,r"""1""",r"""GENCLS""",1, Hshed)
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.strt(0,r"""C:\Users\DELL\Desktop\PhD\energy function\NRPG\energy_function.out""")
psspy.run(0, 1.0,0,1,0)
psspy.dist_bus_fault(157,1, 220.0,[0.0,-0.2E+10])
psspy.run(0, tf,0,1,0)
psspy.dist_clear_fault(1)
psspy.run(0, delay,0,1,0)
psspy.dist_machine_trip(22,r"""2""")
psspy.dist_machine_trip(22,r"""1""")
psspy.run(0, 10.0,0,1,0)
