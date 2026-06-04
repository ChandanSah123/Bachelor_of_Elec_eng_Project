# File:"C:\Users\DELL\Desktop\NRPG_shed1.py", generated on WED, DEC 30 2015  12:52, release 32.00.03
P=360;
Xd=0.2297;
H=3.7;
shed=0.9;
remain=1-shed;
Premain=P*remain
Pshed=P*shed
Xdremain=Xd*(1+shed/remain)
Xdshed=Xd*(1+remain/shed)
Hremain=H*remain
Hshed=H*shed

psspy.case(r"""C:\Users\DELL\Desktop\PhD\Impact of Latency on TSC\NRPG\NRexciter_pristine_cls.sav""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\DELL\Desktop\PhD\Impact of Latency on TSC\NRPG\NRexciter_pristine_cls.dyr""","","","")
psspy.machine_data_2(22,r"""1""",[_i,_i,_i,_i,_i,_i],[ Premain,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_data_2(22,r"""2""",[_i,_i,_i,_i,_i,_i],[ Pshed,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_data_2(22,r"""1""",[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xdremain,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_data_2(22,r"""2""",[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xdshed,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.change_plmod_con(22,r"""1""",r"""GENCLS""",1, Hremain)
psspy.change_plmod_con(22,r"""2""",r"""GENCLS""",1, Hshed)
psspy.save(r"""C:\Users\DELL\Desktop\PhD\Impact of Latency on TSC\NRPG\NRexciter_pristine_cls.sav""")
psspy.dyda(0,1,[2,1,0],0,r"""C:\Users\DELL\Desktop\PhD\Impact of Latency on TSC\NRPG\NRexciter_pristine_cls.dyr""")
psspy.dynamics_solution_param_2([_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f, 0.001,_f,_f,_f,_f,_f])
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.strt(0,r"""C:\Users\DELL\Desktop\NRPG.out""")
psspy.run(0, 1.0,0,1,0)
psspy.dist_bus_fault(157,1, 220.0,[0.0,-0.2E+10])
psspy.run(0, 1.392,0,1,0)
psspy.dist_clear_fault(1)
psspy.run(0, 1.56,0,1,0)
psspy.dist_machine_trip(22,r"""2""")
psspy.run(0, 10.0,0,1,0)
