# File:"C:\Users\DELL\Desktop\nrgrid\busfault.py", generated on TUE, NOV 18 2014  23:08, release 32.00.03
b=157

psspy.case(r"""C:\Users\DELL\Desktop\nrgrid\NRexciter_pristine.sav""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\DELL\Desktop\nrgrid\NRexciter_pristine.dyr""","","","")
psspy.resq(r"""C:\Users\DELL\Desktop\nrgrid\NRexciter_pristine.seq""")
psspy.fdns([0,0,0,1,1,0,99,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.strt(0,r"""C:\Users\DELL\Desktop\nrgrid\allbus.out""")
psspy.run(0, 1.0,100,1,0)
psspy.dist_bus_fault(b,1, 220.0,[0.0,-0.2E+10])
psspy.run(0, 1.38,100,1,0)
psspy.dist_clear_fault(1)
psspy.run(0, 2.0,100,1,0)
