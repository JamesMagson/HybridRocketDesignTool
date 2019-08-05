import tkinter as tk
from tkinter import ttk
import math
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile 


import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
from numpy import sin, cos, pi, tan


#initial estimates
#inputs

#in seconds
total_burn_time = 7
# in m/s
exhaust_velocity = 2288.0844

#earth gravitational parameter
e_g_p = 398600500.00
#earth mean radius
e_m_r = 6378.14

#in meters
launch_altitude = 0
target_altitude = 3048 

#empty mass target in kg
empty_mass_target = 15

#in kg/s
start_flow_rate = 0.53

#ullage estimate
ullage = 0.12

LARGE_FONT = ("Verdana",12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana",8)

def popupmsg(msg):
    popup = tk.Tk()
    
    def leavemini():
        popup.destroy()
    
    popup.wm_title("!")
    label = ttk.Label(popup,text=msg,font=NORM_FONT)
    label.pack(side="top", fill = "x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = leavemini)
    B1.pack()
    popup.mainloop()


def file_save():
    df = pd.DataFrame({'Initial Estimate':[exvel2,lalt2,talt2,tb2,emt2,snfr2,0,0,0,0,0],
                       'Eng Model':[start_thrust2,0,0,0,0,0,0,0,0,0,0],
                       'Injector':[ihd2,of2,tpcp2,lc2,doo2,mfrn2,mpd2,Pr2,Ys2,roi2,SF2],
                       'Combustion Chamber':[avg_rate2,vap_avg_rate2,fuel_density2,max_flux2,0,0,0,0,0,0,0],
                       'Nozzle':[expanr1,tr1,eia1,tn1,te1,0,0,0,0,0,0]})
    writer = ExcelWriter('data.xlsx')
    df.to_excel(writer,'Sheet1',index=False)
    writer.save()
 
def open_file():
    df = pd.read_excel('data.xlsx',sheetname='Sheet1')
    initentry = df['Initial Estimate']
    engentry = df['Eng Model']
    injectentry = df['Injector']
    ccentry = df['Combustion Chamber']
    nozentry = df['Nozzle']
    
    exvel.set(initentry[0])
    lalt.set(initentry[1])
    talt.set(initentry[2])
    tb.set(initentry[3])
    emt.set(initentry[4])
    snfr.set(initentry[5])
    
    start_thrust.set(engentry[0])
    
    ihd.set(injectentry[0])
    of.set(injectentry[1])
    tpcp.set(injectentry[2])
    lc.set(injectentry[3])
    doo.set(injectentry[4])
    mfrn.set(injectentry[5])
    mpd.set(injectentry[6])
    Pr.set(injectentry[7])
    Ys.set(injectentry[8])
    roi.set(injectentry[9])
    SF.set(injectentry[10])
    
    avg_rate.set(ccentry[0])
    vap_avg_rate.set(ccentry[1])
    fuel_density.set(ccentry[2])
    max_flux.set(ccentry[3])
    
    expanr.set(nozentry[0])
    tr.set(nozentry[1])
    eia.set(nozentry[2])
    tn.set(nozentry[3])
    te.set(nozentry[4])


class HDT(tk.Tk):
     def __init__(self,*args,**kwargs):
        tk.Tk.__init__(self,*args,**kwargs)
        tk.Tk.iconbitmap(self,default="raptoricon.ico") 
        tk.Tk.wm_title(self, "Hybrid Rocket Engine Design Tool")
        container = tk.Frame(self)
        container.pack(side="top",fill="both",expand =True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff = 0)
        filemenu.add_command(label="Save", command = lambda: file_save())
        filemenu.add_command(label="Open", command = lambda: open_file())
        filemenu.add_separator()
        filemenu.add_command(label="Create Report", command = lambda: popupmsg('Not supported yet'))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command = lambda: tk.Tk.destroy(self))
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_command(label="Home", command = lambda: self.show_frame(StartPage))
        menubar.add_command(label="Initial Estimates", command = lambda: self.show_frame(InitialEstimates))
        menubar.add_command(label="Engine Model", command = lambda: self.show_frame(EngineModel))
        menubar.add_command(label="Injector", command = lambda: self.show_frame(Injector))
        menubar.add_command(label="Combustion Chamber", command = lambda: self.show_frame(CombustionChamber))
        menubar.add_command(label="Nozzle", command = lambda: self.show_frame(Nozzle))
        menubar.add_command(label="Drag Coefficient", command = lambda: self.show_frame(DragCoeff))
        menubar.add_command(label="Atmosphere", command = lambda: self.show_frame(Atmosphere))
        menubar.add_command(label="Simulation", command = lambda: self.show_frame(Simulation))

        tk.Tk.config(self, menu=menubar)
        
        self.frames = {}
        for F in (StartPage, InitialEstimates,EngineModel ,Injector, CombustionChamber, 
                  Nozzle,DragCoeff,Atmosphere,Simulation):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column = 0, sticky="nsew")
        self.show_frame(StartPage)
     def show_frame(self,cont):
        frame = self.frames[cont]
        frame.tkraise()

def qf(param):
    print(param)
    
def initial_estimates(a,a0,a1,a2,emt1,snfr1):
    global a12, a13, snfr2, exvel2,lalt2,talt2,tb2,emt2, a8
    snfr2 = snfr1
    exvel2 = a
    lalt2 = a1
    talt2 = a2
    tb2 = a0
    emt2 = emt1
    a3 = (e_g_p/(e_m_r +(a1/1000))-(e_g_p/(e_m_r+(a2/1000))))
    a4 = (2*a3*1000)**0.5
    a5 = 9.81*a0
    a6 = a5
    a7 = a4+ a5+ a6
    a8 = math.exp(a7/a) if a else 0
    a9 = emt1*(a8-1)
    a10 = snfr1*0.7
    a11 = emt1 + a9
    a12 = a9*(1-0.12)
    a13 = 0.12*a9
    a14 = (2*a12)/(snfr1+a10) if a else 0
    return a3,a4,a5,a6,a7,a8,a9,a11,a12,a13,a14

def injectorcalc(ihd1,of1,tpcp1,lc1,doo1,mfrn1):
    global ihd2,of2,tpcp2,lc2,doo2,mfrn2
    ihd2 = ihd1
    of2 = of1
    tpcp2=tpcp1
    lc2=lc1
    doo2=doo1
    mfrn2=mfrn1
    aoo = (pi*(ihd1)**2)/4.0
    mfrn2 = of1*mfrn1
    noo = mfrn2/((aoo)*(((2*doo1*tpcp1)/lc1)**(0.5))) if aoo else 0
    return aoo,noo

def injectorcalcstress(mpd1,Pr1,Ys1,roi1,SF1):
    global t, mpd2,Pr2,Ys2,roi2,SF2
    mpd2=mpd1
    Pr2=Pr1
    Ys2=Ys1
    roi2=roi1
    SF2=SF1
    t = (((0.375*mpd1*(roi1**3)*(1+Pr1))/((Ys1)/SF1))**0.5)*1000 if Ys1 else 0
    return t

def engcalc(start_thrust1):
    global start_thrust2,liquid_burn_time,vapour_burn_time
    start_thrust2 = start_thrust1
    burnout_thrust = 0.7* start_thrust1
    burnout_flow_rate = snfr2*0.7
    liquid_burn_time = (2*a12)/(snfr2+burnout_flow_rate)
    liquid_impulse = (start_thrust1 + burnout_thrust)*0.5*liquid_burn_time
    flow_grad = (burnout_flow_rate - snfr2)/liquid_burn_time
    thrust_grad = (burnout_thrust- start_thrust1)/liquid_burn_time
    vapour_impulse = liquid_impulse*0.15 
    vapour_burn_time = (2*vapour_impulse)/burnout_thrust 
    total_burn_time = liquid_burn_time + vapour_burn_time 
    total_impulse = liquid_impulse + vapour_impulse 
    
    result = (start_thrust1,burnout_thrust, exvel2, snfr2, burnout_flow_rate, a12,a13,liquid_burn_time
              ,liquid_impulse, flow_grad, thrust_grad,vapour_impulse,vapour_burn_time,total_burn_time,total_impulse)
    
    return result

def portgeom(avg_rate1, vap_avg_rate1, fuel_density1, max_flux1):
    global avg_rate2,vap_avg_rate2,fuel_density2,max_flux2
    avg_rate2=avg_rate1
    vap_avg_rate2=vap_avg_rate1
    fuel_density2=fuel_density1
    max_flux2=max_flux1
    nitrous_mass = a12*(7.0/8.0)
    fuel_mass = a12*(1.0/8.0)
    end_flow_rate = 0.7*snfr2
    init_port_area = snfr2/max_flux1
    init_port_diameter = (((4*init_port_area)/pi)**0.5)*1000
    final_liquid_port_diameter = init_port_diameter + ((liquid_burn_time*avg_rate1)*2)
    final_port_diameter = final_liquid_port_diameter + ((vapour_burn_time*vap_avg_rate1)*2)
    mid_burn_diameter = (init_port_diameter+final_port_diameter)/2
    mid_flow_rate = (snfr2+(0.7*snfr2))/2
    length_grain = mid_flow_rate/(8*(avg_rate1/1000)*fuel_density1*pi*(mid_burn_diameter/1000))
    port_length = length_grain + final_liquid_port_diameter/1000
    result = (nitrous_mass,fuel_mass,end_flow_rate,init_port_area,init_port_diameter
              ,final_liquid_port_diameter,final_port_diameter,mid_burn_diameter,mid_flow_rate,length_grain
              ,port_length)
    return result


#def sim(fuselage_diameter1):
    
    
    
    
    



class StartPage(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        label=tk.Label(self,text="Home", font = LARGE_FONT)
        label.grid(row=0,column=0)
        
        
class InitialEstimates(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
#       #object is label, bg background, fg foreground
#       #input labels
        ilabel1 = tk.Label(self, text="Exhaust Velocity (m/s)", fg="black").grid(row=4,sticky = "e")
        ilabel2 = tk.Label(self, text="Launch Alitutde (m)", fg="black").grid(row=5,sticky = "e")
        ilabel3 = tk.Label(self, text="Target Alitutde (m)", fg="black").grid(row=6,sticky = "e")
        ilabel4 = tk.Label(self, text="Estimate burn time (s)", fg="black").grid(row=7,sticky = "e")
        ilabel5 = tk.Label(self, text="Empty mass target (kg)", fg="black").grid(row=8,sticky = "e")
        ilabel6 = tk.Label(self, text="Start nozzle flow rate (kg/s)", fg="black").grid(row=9,sticky = "e")
        
#       #output labels 
        olabel1 = tk.Label(self, text="Specific gravitational potential (kJ)").grid(row=11,sticky = "e")
        olabel2 = tk.Label(self, text="Velocity (m/s)").grid(row=12,sticky = "e")
        olabel3 = tk.Label(self, text="Gravity loss (m/s)").grid(row=13,sticky = "e")
        olabel4 = tk.Label(self, text="Drag loss (m/s)").grid(row=14,sticky = "e")
        olabel5 = tk.Label(self, text="Total velocity (m/s)").grid(row=15,sticky = "e")
        olabel6 = tk.Label(self, text="Mass ratio").grid(row=16,sticky = "e")
        olabel7 = tk.Label(self, text="Propellant mass (kg)").grid(row=17,sticky = "e")
        olabel8 = tk.Label(self, text="Total mass (kg)").grid(row=18,sticky = "e")
        olabel9 = tk.Label(self, text="Start propellant mass (kg)").grid(row=19,sticky = "e")
        olabel10 = tk.Label(self, text="Final vapour mass (kg)").grid(row=20,sticky = "e")
        olabel11 = tk.Label(self, text="Liquid burn time (s)").grid(row=21,sticky = "e")
        global exvel
        global lalt
        global talt
        global tb
        global emt
        global snfr
        exvel = tk.DoubleVar()
        lalt = tk.DoubleVar()
        talt = tk.DoubleVar()
        tb = tk.DoubleVar()
        emt = tk.DoubleVar()
        snfr = tk.DoubleVar()
        
        self.entry1 = tk.Entry(self, textvariable = exvel)
        self.entry1.grid(row=4,column=1)
        self.entry2 = tk.Entry(self, textvariable = lalt)
        self.entry2.grid(row=5,column=1)
        self.entry3 = tk.Entry(self, textvariable = talt)
        self.entry3.grid(row=6,column=1)
        self.entry4 = tk.Entry(self, textvariable = tb)
        self.entry4.grid(row=7,column=1)
        self.entry5 = tk.Entry(self, textvariable = emt)
        self.entry5.grid(row=8,column=1)
        self.entry6 = tk.Entry(self, textvariable = snfr)
        self.entry6.grid(row=9,column=1)
        button1 = ttk.Button(self,text="Calculate",command=self.initial_estimates1).grid(row=10,column=0)
    def initial_estimates1(self):
        a = float(self.entry1.get())
        a0 = float(self.entry4.get())
        a1 = float(self.entry2.get())
        a2 = float(self.entry3.get())
        emt1 = float(self.entry5.get())
        snfr1 = float(self.entry6.get())
        result = initial_estimates(a,a0,a1,a2,emt1,snfr1)
        #print("%d a = %d 2" % (a0,result))
        self.label4 = tk.Label(self, text ="%.2f" % result[0])
        self.label4.grid(row = 11, column=1)
        self.label5 = tk.Label(self, text ="%.2f" % result[1])
        self.label5.grid(row = 12, column =1)
        self.label6 = tk.Label(self, text ="%.2f" % result[2])
        self.label6.grid(row=13,column=1)
        self.label7 = tk.Label(self, text ="%.2f" % result[3])
        self.label7.grid(row = 14, column =1)
        self.label8 = tk.Label(self, text ="%.2f" % result[4])
        self.label8.grid(row = 15, column =1)
        self.label9 = tk.Label(self, text ="%.2f" % result[5])
        self.label9.grid(row = 16, column =1)
        self.label10 = tk.Label(self, text ="%.2f" % result[6])
        self.label10.grid(row = 17, column =1)
        self.label11 = tk.Label(self, text ="%.2f" % result[7])
        self.label11.grid(row = 18, column =1)
        self.label12 = tk.Label(self, text ="%.2f" % result[8])
        self.label12.grid(row = 19, column =1)
        self.label13 = tk.Label(self, text ="%.2f" % result[9])
        self.label13.grid(row = 20, column =1)
        self.label13 = tk.Label(self, text ="%.2f" % result[10])
        self.label13.grid(row = 21, column =1)
    


class EngineModel(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        
        ilabel1 = tk.Label(self, text="Start Thrust (N)", fg="black").grid(row=1,sticky = "e")
        
        olabel1 = tk.Label(self, text="Start Thrust (N)").grid(row=3,column=0)
        olabel2 = tk.Label(self, text="Burnout Thrust (N)").grid(row=3,column=1)
        olabel3 = tk.Label(self, text="Exhaust Velocity(m/s)").grid(row=3,column=2)
        olabel4 = tk.Label(self, text="Start Mass Flow Rate (kg/s)").grid(row=3,column=3)
        olabel5 = tk.Label(self, text="Burnout Mass Flow Rate (kg/s)").grid(row=3,column=4)
        olabel6 = tk.Label(self, text="Start Propellant Mass (kg)").grid(row=3,column=5)
        olabel7 = tk.Label(self, text="Final Vapour Mass(kg)").grid(row=5,column=0)
        olabel8 = tk.Label(self, text="Liquid Burn Time(s)").grid(row=5,column=1)
        olabel9 = tk.Label(self, text="Liquid Impulse (Ns)").grid(row=5,column=2)
        olabel10 = tk.Label(self, text="Mass Gradient").grid(row=5,column=3)
        olabel11 = tk.Label(self, text="Thrust Gradient").grid(row=5,column=4)
        olabel12 = tk.Label(self, text="Vapour Impulse (Ns)").grid(row=5,column=5)
        olabel13 = tk.Label(self, text="Vapour Burn Time (s)").grid(row=5,column=6)
        olabel14 = tk.Label(self, text="Total Burn Time (s)").grid(row=7,column=0)
        olabel15 = tk.Label(self, text="Total Impulse (Ns)").grid(row=7,column=1)
        
        global start_thrust
        start_thrust = tk.DoubleVar()
        
        self.entry1 = tk.Entry(self, textvariable = start_thrust)
        self.entry1.grid(row=1,column=1)
        
        
        button1 = ttk.Button(self,text="Calculate",command= self.engcalc1).grid(row=2,column=0)
    def engcalc1(self):
        start_thrust1 = float(self.entry1.get())
        result = engcalc(start_thrust1)
        self.label1 = tk.Label(self, text ="%.2f" % result[0])
        self.label1.grid(row = 4, column=0)
        self.label2 = tk.Label(self, text ="%.2f" % result[1])
        self.label2.grid(row = 4, column =1)
        self.label3 = tk.Label(self, text ="%.2f" % result[2])
        self.label3.grid(row = 4, column =2)
        self.label4 = tk.Label(self, text ="%.2f" % result[3])
        self.label4.grid(row = 4, column =3)
        self.label5 = tk.Label(self, text ="%.2f" % result[4])
        self.label5.grid(row = 4, column=4)
        self.label6 = tk.Label(self, text ="%.2f" % result[5])
        self.label6.grid(row = 4, column=5)
        self.label7 = tk.Label(self, text ="%.2f" % result[6])
        self.label7.grid(row = 6, column=0)
        self.label8 = tk.Label(self, text ="%.2f" % result[7])
        self.label8.grid(row = 6, column=1)
        self.label9 = tk.Label(self, text ="%.2f" % result[8])
        self.label9.grid(row = 6, column=2)
        self.label10 = tk.Label(self, text ="%.2f" % result[9])
        self.label10.grid(row = 6, column=3)
        self.label11 = tk.Label(self, text ="%.2f" % result[10])
        self.label11.grid(row = 6, column=4)
        self.label12 = tk.Label(self, text ="%.2f" % result[11])
        self.label12.grid(row = 6, column=5)
        self.label13 = tk.Label(self, text ="%.2f" % result[12])
        self.label13.grid(row = 6, column=6)
        self.label14 = tk.Label(self, text ="%.2f" % result[13])
        self.label14.grid(row = 8, column=0)
        self.label15 = tk.Label(self, text ="%.2f" % result[14])
        self.label15.grid(row = 8, column=1)
        
        
        time_array = np.linspace(0,result[7],100)
        thrust_array = np.array([])
        t0 = 0
        for t in time_array:
            t0 = result[0]+t*result[10]
            thrust_array =np.append(thrust_array,t0)
        
        vapour_time = np.array([time_array[-1],result[13]])
        vapour_thrust = np.array([thrust_array[-1],0])
        
        flow_array = np.array([])
        f0 = 0
        for f in time_array:
            f0 = result[3]+f*result[9]
            flow_array =np.append(flow_array,f0)
        
        vapour_flow = np.array([flow_array[-1],0])
        print(vapour_flow)
        
        
        f1 = Figure(figsize=(5,5), dpi=100)
       
        
        canvas = FigureCanvasTkAgg(f1,self)
        canvas.get_tk_widget().grid(row=9,column=1,columnspan=3,rowspan=20)
#        canvas._tkcanvas.grid()
        a = f1.add_subplot(111)
        a.plot(time_array,thrust_array,label='Liquid Phase')
        a.plot(vapour_time,vapour_thrust, label='Vapour Phase')
        a.set_xlabel('Time (s)')
        a.set_ylabel('Thrust (N)')
        a.legend(loc='upper right')
        a.grid(b=True, which='major',color='black',linewidth=0.2)
        canvas.draw()
        
        toolbarFrame = tk.Frame(self)
        toolbarFrame.grid(row=31,column=1)
        toolbar = NavigationToolbar2TkAgg(canvas,toolbarFrame)
        toolbar.update()
  
#        f2 = Figure(figsize=(5,5), dpi=100)
#       
#        
#        canvas = FigureCanvasTkAgg(f2,self)
#        canvas.get_tk_widget().grid(row=9,column=3,columnspan=3,rowspan=20)
#        canvas._tkcanvas.grid()
#        a = f1.add_subplot(111)
#        a.plot(time_array,flow_array,vapour_time,vapour_flow)
#        a.set_xlabel('Time (s)')
#        a.set_ylabel('Flow Rate (kg/s)')
#        canvas.show()
#        
#        toolbarFrame = tk.Frame(self)
#        toolbarFrame.grid(row=30,column=0)
#        toolbar = NavigationToolbar2TkAgg(canvas,toolbarFrame)
#        toolbar.update()
#        
class Injector(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        #geometry
        ilabel1 = tk.Label(self, text="Geometry", fg="black").grid(row=4,sticky = "e")
        ilabel2 = tk.Label(self, text="Injector hole diameter (mm)", fg="black").grid(row=5,sticky = "e")
        ilabel3 = tk.Label(self, text="Oxidiser/Fuel ratio" ,fg="black",).grid(row=6,sticky = "e")
        ilabel4 = tk.Label(self, text="Tank pressure - Chamber pressure (Bar)", fg="black").grid(row=7,sticky = "e")
        ilabel5 = tk.Label(self, text="Loss coefficient", fg="black").grid(row=8,sticky = "e")
        ilabel6 = tk.Label(self, text="Density of oxidiser (kg/m^3)", fg="black").grid(row=9,sticky = "e")
        ilabel7 = tk.Label(self, text="Mass flow rate nitrous (kg/s)", fg="black").grid(row=10,sticky = "e")
        
        
        #stress
        ilabel9 = tk.Label(self, text="Stress Analysis", fg="black").grid(row=12,sticky = "e")
        ilabel10 = tk.Label(self, text="Max pressure drop (Pa)", fg="black").grid(row=13,sticky = "e")
        ilabel11 = tk.Label(self, text="Poisson Ratio", fg="black").grid(row=14,sticky = "e")
        ilabel12 = tk.Label(self, text="Yield Stress (Pa)", fg="black").grid(row=15,sticky = "e")
        ilabel13 = tk.Label(self, text="Radius of injector (m)", fg="black").grid(row=16,sticky = "e")
        ilabel14 = tk.Label(self, text="Safety Factor", fg="black").grid(row=17,sticky = "e")
        
        global ihd
        global of
        global tpcp
        global lc
        global doo
        global mfrn
        global mpd
        global Pr
        global Ys
        global roi 
        global SF
        
        ihd = tk.DoubleVar()
        of = tk.DoubleVar()
        tpcp = tk.DoubleVar()
        lc = tk.DoubleVar()
        doo = tk.DoubleVar()
        mfrn = tk.DoubleVar()
        mpd = tk.DoubleVar()
        Pr = tk.DoubleVar()
        Ys = tk.DoubleVar()
        roi = tk.DoubleVar()
        SF = tk.DoubleVar()
        
        self.entry1 = tk.Entry(self, textvariable = ihd)
        self.entry1.grid(row=5,column=1)
        self.entry2 = tk.Entry(self, textvariable = of)
        self.entry2.grid(row=6,column=1)
        self.entry3 = tk.Entry(self, textvariable = tpcp)
        self.entry3.grid(row=7,column=1)
        self.entry4 = tk.Entry(self, textvariable = lc)
        self.entry4.grid(row=8,column=1)
        self.entry5 = tk.Entry(self, textvariable = doo)
        self.entry5.grid(row=9,column=1)
        self.entry6 = tk.Entry(self, textvariable = mfrn)
        self.entry6.grid(row=10,column=1)
        self.entry7 = tk.Entry(self, textvariable = mpd)
        self.entry7.grid(row=13,column=1)
        self.entry8 = tk.Entry(self, textvariable = Pr)
        self.entry8.grid(row=14,column=1)
        self.entry9 = tk.Entry(self, textvariable = Ys)
        self.entry9.grid(row=15,column=1)
        self.entry10 = tk.Entry(self, textvariable = roi)
        self.entry10.grid(row=16,column=1)
        self.entry11 = tk.Entry(self, textvariable = SF)
        self.entry11.grid(row=17,column=1)
        
        olabel1 = tk.Label(self, text="Area of orifice (m^2)").grid(row=5,column=2,sticky = "e")
        olabel2 = tk.Label(self, text="Number of orificies").grid(row=6,column=2,sticky = "e")
        olabel3 = tk.Label(self, text="Thickness (mm)").grid(row=13,column=2,sticky = "e")


        button1 = tk.Button(self,text="Calculate",command=self.injectorcalc1).grid(row=11,column=0)
        button2 = tk.Button(self,text="Calculate",command=self.injectorcalcstress1).grid(row=18,column=0)
        
    def injectorcalc1(self):
        ihd1 = float(self.entry1.get())
        of1 = float(self.entry2.get())
        tpcp1 = float(self.entry3.get())
        lc1 = float(self.entry4.get())
        doo1 = float(self.entry5.get())
        mfrn1 = float(self.entry6.get())
        result = injectorcalc(ihd1,of1,tpcp1,lc1,doo1,mfrn1)
        #print("%d a = %d 2" % (a0,result))
        self.label4 = tk.Label(self, text ="%.2f" % result[0])
        self.label4.grid(row = 5, column=3)
        self.label5 = tk.Label(self, text ="%.2f" % result[1])
        self.label5.grid(row = 6, column =3)
        
    def injectorcalcstress1(self):
        mpd1 = float(self.entry7.get())
        Pr1 = float(self.entry8.get())
        Ys1 = float(self.entry9.get())
        roi1 = float(self.entry10.get())
        SF1 = float(self.entry11.get())
        result = injectorcalcstress(mpd1,Pr1,Ys1,roi1,SF1)
        self.label4 = tk.Label(self, text ="%.2f" % result)
        self.label4.grid(row = 13, column=3)
        
        
class CombustionChamber(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
#        label1=tk.Label(self,text="Combustion Chamber", font = LARGE_FONT)
#        label1.grid(row=0,column=0)
        
        ilabel1 = tk.Label(self, text="Average regression rate (mm/s)", fg="black").grid(row=0,sticky = "e")
        ilabel2 = tk.Label(self, text="Vapour average regression rate (mm/s)", fg="black").grid(row=1,sticky = "e")
        ilabel3 = tk.Label(self, text="Density of fuel (kg/m^3)" ,fg="black",).grid(row=2,sticky = "e")  
        ilabel4 = tk.Label(self, text="Maximum Flux (kg/s/m^2)" ,fg="black",).grid(row=3,sticky = "e")  
        
        global avg_rate
        global vap_avg_rate
        global fuel_density
        global max_flux
        
        avg_rate = tk.DoubleVar()
        vap_avg_rate = tk.DoubleVar()
        fuel_density = tk.DoubleVar()
        max_flux = tk.DoubleVar()
        
             
        self.entry1 = tk.Entry(self, textvariable = avg_rate)
        self.entry1.grid(row=0,column=1)
        self.entry2 = tk.Entry(self, textvariable = vap_avg_rate)
        self.entry2.grid(row=1,column=1)
        self.entry3 = tk.Entry(self, textvariable = fuel_density)
        self.entry3.grid(row=2,column=1)
        self.entry4 = tk.Entry(self, textvariable = max_flux)
        self.entry4.grid(row=3,column=1)
        
        button1 = tk.Button(self,text="Calculate",command=self.portgeom1).grid(row=4,column=0)

        olabel1 = tk.Label(self, text="Oxidiser mass (kg)").grid(row=5,column=0,sticky = "e")
        olabel2 = tk.Label(self, text="Fuel mass (kg)").grid(row=6,column=0,sticky = "e")
        olabel3 = tk.Label(self, text="End flow rate (kg/s)").grid(row=7,column=0,sticky = "e")
        olabel4 = tk.Label(self, text="Initial port area (m^2)").grid(row=8,column=0,sticky = "e")
        olabel5 = tk.Label(self, text="Initial port diameter (m)").grid(row=9,column=0,sticky = "e")
        olabel6 = tk.Label(self, text="Final liquid port diameter(mm)").grid(row=10,column=0,sticky = "e")
        olabel7 = tk.Label(self, text="Final port diameter(mm)").grid(row=11,column=0,sticky = "e")
        olabel8 = tk.Label(self, text="Mid burn port diamater (mm)").grid(row=12,column=0,sticky = "e")
        olabel9 = tk.Label(self, text="Mid flow rate (kg/s)").grid(row=13,column=0,sticky = "e")
        olabel10 = tk.Label(self, text="Grain length (m)").grid(row=14,column=0,sticky = "e")
        olabel11 = tk.Label(self, text="Port length (m)").grid(row=15,column=0,sticky = "e")
        
    def portgeom1(self):
        avg_rate1 = float(self.entry1.get())
        vap_avg_rate1 = float(self.entry2.get())
        fuel_density1 = float(self.entry3.get())
        max_flux1 = float(self.entry4.get())
        result = portgeom(avg_rate1,vap_avg_rate1,fuel_density1,max_flux1)
        
        self.label1 = tk.Label(self, text ="%.2f" % result[0])
        self.label1.grid(row = 5, column=1)
        self.label2 = tk.Label(self, text ="%.2f" % result[1])
        self.label2.grid(row = 6, column=1)
        self.label3 = tk.Label(self, text ="%.2f" % result[2])
        self.label3.grid(row = 7, column=1)
        self.label4 = tk.Label(self, text ="%.4f" % result[3])
        self.label4.grid(row = 8, column=1)
        self.label5 = tk.Label(self, text ="%.2f" % result[4])
        self.label5.grid(row = 9, column=1)
        self.label6 = tk.Label(self, text ="%.2f" % result[5])
        self.label6.grid(row = 10, column=1)
        self.label7 = tk.Label(self, text ="%.2f" % result[6])
        self.label7.grid(row = 11, column=1)
        self.label8 = tk.Label(self, text ="%.2f" % result[7])
        self.label8.grid(row = 12, column=1)
        self.label9 = tk.Label(self, text ="%.2f" % result[8])
        self.label9.grid(row = 13, column=1)
        self.label10 = tk.Label(self, text ="%.2f" % result[9])
        self.label10.grid(row = 14, column=1)
        self.label11 = tk.Label(self, text ="%.2f" % result[10])
        self.label11.grid(row = 15, column=1)

class Nozzle(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        #label1=tk.Label(self,text="Nozzle", font = LARGE_FONT)
        #label1.grid(row=0,column=0)
        #ilabel1 = tk.Label(self, text="Exhaust Velocity (m/s)", fg="black").pack(side = tk.LEFT)
        #ilabel1 = tk.Label(self, text="Exhaust Velocity (m/s)", fg="black").pack(side = tk.LEFT)
        
        ilabel1 = tk.Label(self, text="Expansion Ratio", fg="black").grid(row=0,sticky = "e")
        ilabel2 = tk.Label(self, text="Throat Radius (mm)", fg="black").grid(row=1,sticky = "e")
        ilabel3 = tk.Label(self, text="Entrant Initial Angle (rad)" ,fg="black",).grid(row=2,sticky = "e")
        ilabel4 = tk.Label(self, text="ThetaN", fg="black").grid(row=3,sticky = "e")
        ilabel5 = tk.Label(self, text="ThetaE", fg="black").grid(row=4,sticky = "e")
        
        global expanr
        global tr
        global eia
        global tn
        global te 
        
        expanr = tk.DoubleVar()
        tr = tk.DoubleVar()
        eia = tk.DoubleVar()
        tn = tk.DoubleVar()
        te = tk.DoubleVar()
        
        self.entry1 = tk.Entry(self, textvariable = expanr)
        self.entry1.grid(row=0,column=1)
        self.entry2 = tk.Entry(self, textvariable = tr)
        self.entry2.grid(row=1,column=1)
        self.entry3 = tk.Entry(self, textvariable = eia)
        self.entry3.grid(row=2,column=1)
        self.entry4 = tk.Entry(self, textvariable = tn)
        self.entry4.grid(row=3,column=1)
        self.entry5 = tk.Entry(self, textvariable = te)
        self.entry5.grid(row=4,column=1)
        
        
        button1 = tk.Button(self,text="Calculate",command=self.nozzle_calc).grid(row=5,column=0)
    def nozzle_calc(self):
    
        e = float(self.entry1.get())
        Rt = float(self.entry2.get())
        entrant_initial_angle = float(self.entry3.get())
        theta_N = float(self.entry4.get())
        theta_E = float(self.entry5.get())
        
        global expanr1,tr1,eia1,tn1,te1
        
        expanr1 = e
        tr1 = Rt
        eia1 = entrant_initial_angle
        tn1 = theta_N
        te1 = theta_E
    
        
        #entrant functions
        entrant_angles = np.linspace(entrant_initial_angle, -pi/2,100)
        
        xe = np.array([])
        xe1 = 0
        for i in entrant_angles:
            xe1 = 1.5*Rt*np.cos(i)
            xe = np.append(xe,xe1)
        
        ye = np.array([])
        ye1 = 0
        for j in entrant_angles:
            ye1 = 1.5*Rt*np.sin(j) + 2.5*Rt
            ye = np.append(ye,ye1)
        
        #exit section
        exit2 = theta_N - pi/2
        exit_angles = np.linspace(-pi/2, exit2,100)
        
        xexit = np.array([])
        xexit1 = 0
        for i in exit_angles:
            xexit1 = 0.382*Rt*cos(i)
            xexit = np.append(xexit, xexit1)
            
        yexit = np.array([])
        yexit1 = 0
        for j in exit_angles:
            yexit1 = 0.382*Rt*sin(j) + 1.382*Rt
            yexit = np.append(yexit,yexit1)
        
        #bell section
        Nx = 0.382*Rt*cos(theta_N-pi/2)
        Ny = 0.382*Rt*sin(theta_N-pi/2) + 1.382*Rt
        Ex = 0.8*(((e**0.5)-1)*Rt)/(tan(15*pi/180))
        Ey = ((e)**0.5)*Rt
        m1 = tan(theta_N)
        m2 = tan(theta_E)
        C1 = Ny - m1*Nx
        C2 = Ey - m2*Ex 
        Qx = (C2 - C1)/(m1 - m2)
        Qy = (m1*C2 - m2*C1)/(m1 - m2)
        t = np.linspace(0,1,100)
        
        xbell = np.array([])
        xbell1 = 0
        for i in t:
            xbell1 = ((1-i)**2)*Nx + 2*(1-i)*i*Qx + (i**2)*Ex
            xbell = np.append(xbell, xbell1)
        
        ybell = np.array([])
        ybell1 = 0
        for j in t:
            ybell1 = ((1-j)**2)*Ny + 2*(1-j)*j*Qy + (j**2)*Ey
            ybell = np.append(ybell, ybell1)
        #print(xe,ye)
        f = Figure(figsize=(5,5), dpi=100)
       
        
        canvas = FigureCanvasTkAgg(f,self)
        canvas.get_tk_widget().grid(row=6,column=0,columnspan=3,rowspan=20)
#        canvas._tkcanvas.grid()
        a = f.add_subplot(111)
        a.plot(xe,ye,xexit,yexit,xbell,ybell)
        a.set_xlabel('X (mm)')
        a.set_ylabel('Y (mm)')
        a.grid(b=True, which='major',color='black',linewidth=0.2)
        canvas.draw()
        
        toolbarFrame = tk.Frame(self)
        toolbarFrame.grid(row=27,column=0)
        toolbar = NavigationToolbar2TkAgg(canvas,toolbarFrame)
        toolbar.update()
        
class DragCoeff(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        
        mach_array = [0.0,0.2,0.4,0.6,0.70,0.75,0.8,0.85,0.9,
                      0.95,1.0,1.07,1.1,1.15,1.2,1.3,1.5,1.75,
                      2.2,2.5,3.0,3.5,4.0,4.5,5.0,8.0]
        poweron_array = [0.476,0.458,0.444,0.435,0.434,0.434,0.434,0.442,0.460,0.488,
                         0.516,0.551,0.556,0.556,0.550,0.535,0.51,0.488,0.452,0.433,0.407,
                         0.388,0.374,0.365,0.360,0.350]
        poweroff_array = [0.548,0.518,0.503,0.496,0.498,0.500,0.508,0.517,0.536,0.570,
                          0.632,0.727,0.743,0.733,0.718,0.693,0.661,0.624,0.563,0.528,
                          0.475,0.432,0.405,0.385,0.375,0.350]
        
        f = Figure(figsize=(5,5), dpi=100)
       
        
        canvas = FigureCanvasTkAgg(f,self)
        canvas.get_tk_widget().grid(row=1,column=1,columnspan=3,rowspan=20)

        a = f.add_subplot(111)
        a.plot(mach_array,poweron_array, label='Power On')
        a.plot(mach_array,poweroff_array, label = 'Power Off')
        a.set_xlabel('Mach Number')
        a.set_ylabel('Drag Coefficient')
        a.legend(loc='upper right')
        a.grid(b=True, which='major',color='black',linewidth=0.2)
        canvas.draw()
        
        toolbarFrame = tk.Frame(self)
        toolbarFrame.grid(row=28,column=1)
        toolbar = NavigationToolbar2TkAgg(canvas,toolbarFrame)
        toolbar.update()
        

        

class Atmosphere(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        global height_array
        global atm_density_array
        global sound_speed_array
        global atm_density_log_array 
        
        height_array =[0,50,100,200,300,400,500,1000,1500,2000,3000,4000,5000,7000,
        10000,15000,20000,25000,30000,35000,40000,45000,50000,60000,
        70000,80000,90000,100000,120000,150000,200000,300000,500000,
        1000000]
        atm_density_array = [1.225E+00,1.219E+00,1.213E+00,1.202E+00,1.190E+00,1.179E+00,
        1.167E+00,1.112E+00,1.058E+00,1.007E+00,9.093E-01,8.194E-01,
        7.365E-01,5.900E-01,4.135E-01,1.948E-01,8.891E-02,4.008E-02,
        1.841E-02,8.463E-03,3.996E-03,1.966E-03,1.027E-03,3.097E-04,
        8.282E-05,1.846E-05,3.416E-06,5.607E-07,2.222E-08,2.076E-09,
        2.571E-10,1.916E-11,5.215E-13,3.560E-15]
        atm_density_log_array = [0.2029,0.1981,0.1933,0.1837,0.1740,0.1644,0.1547,0.1059,0.0565,
        0.0066,-0.0951,-0.1992,-0.3059,-0.5276,-0.8831,-1.6360,-2.4201,
        -3.2168,-3.9949,-4.7721,-5.5225,-6.2318,-6.8811,-8.0799,-9.3988,
        -10.8999,-12.5870,-14.3941,-17.6223,-19.9928,-22.0816,-24.6781,
        -28.2821,-33.2691]
        sound_speed_array = [340.30,340.10,339.91,339.53,339.14,338.76,338.40,336.40,334.49,
        332.50,328.60,324.60,320.60,312.30,299.50,295.10,295.10,298.40,
        301.70,308.30,317.20,325.80,329.80,315.10,297.10,282.50,1000000,
        1000000,1000000,1000000,1000000,1000000,1000000,1000000]
        
        f = Figure(figsize=(5,5), dpi=100)
        canvas = FigureCanvasTkAgg(f,self)
        canvas.get_tk_widget().grid(row=1,column=1,columnspan=3,rowspan=20)
        a = f.add_subplot(111)
        a.plot(height_array[0:23],atm_density_array[0:23])
        a.set_xlabel('Height (m)')
        a.set_ylabel('Atmospheric Density (kg/m^3)')
        a.grid(b=True, which='major',color='black',linewidth=0.2)
        canvas.draw()
        

class Simulation(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
                
        ilabel1 = tk.Label(self, text="Fuselage Diameter (m)", fg="black").grid(row=0,sticky = "e")
        
        global fuselage_diameter

        fuselage_diameter = tk.DoubleVar()        
        
        self.entry1 = tk.Entry(self, textvariable = fuselage_diameter)
        self.entry1.grid(row=0,column=1)
        
        button1 = tk.Button(self,text="Calculate",command=lambda: popupmsg('Not supported yet')).grid(row=5,column=0)
    
    def sim1(self):
        fuselage_diameter1 = float(self.entry1.get())
        print(fuselage_diameter1)
      


app = HDT()
app.geometry("980x720")
app.mainloop()
    