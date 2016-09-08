# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:54:32 2016

@author: viherbos
"""

from PIL import Image, ImageTk
from Tkinter import Tk, BOTH, W, N, E, S
from Tkinter import Checkbutton, BooleanVar, IntVar, DoubleVar, Spinbox
from Tkinter import StringVar, Toplevel, Message
from ttk import Frame, Button, Label, Entry
import numpy as np
import matplotlib.pyplot as plt

import find_coeff as fc


class gui(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.parent = parent
        self.initUI()
        self.centerUI(w=600,h=250)

    def initUI(self):
        
        self.parent.title("FIND BLR COEFFICIENT VALUE")
        self.pack(fill=BOTH, expand=True)
        
        self.columnconfigure(0, weight=1)
        # self.rowconfigure(0, weight=1)
        # weight attibute is used to make them growable
        
        ###### GUI Control Variables ######
        self.LIMIT_L            = IntVar()
        self.LIMIT_H            = IntVar()
        self.PULSE_R            = IntVar()
        self.PULSE_L            = IntVar()
        self.pulse_height       = DoubleVar()
        self.hdf5_file          = StringVar()
        self.PMT                = IntVar()
        self.EVENT              = IntVar()
        self.amplitude_range    = DoubleVar()
        self.delta              = DoubleVar()
        self.noise_sigma        = DoubleVar()
        self.coeff              = DoubleVar()
        #self.DRAW               = BooleanVar()
        
        
        search = Image.open("next_logo.jpg")
        search_temp = search.resize((170, 200), Image.ANTIALIAS)
        search_aux = ImageTk.PhotoImage(search_temp)
        label1 = Label(self, image=search_aux)
        label1.image = search_aux
        label1.grid(row=0, column=0,
				columnspan=10, rowspan=10, sticky=E+W+S+N)
        
        
        self.hdf5_file.set("2052.h5.z")
        e1 = Entry(self, textvariable=self.hdf5_file, width=30)
        e1.grid(row=1,column=1, sticky=W, columnspan=5, pady=5)
        e1_label = Label(self, text="HDF5 file")
        e1_label.grid(row=0,column=1,sticky=W, columnspan=5, pady=5)        
               
        self.PMT.set("0")
        sb1 = Spinbox(self, from_=0, to=12, 
                      width=3, textvariable=self.PMT)
        sb1.grid(row=3,column=2, sticky=W)
        sb1_label = Label(self, text="PMT")
        sb1_label.grid(row=2,column=2, padx=0, sticky=W)
        
        self.EVENT.set("0")
        sb1 = Spinbox(self, from_=0, to=1000, 
                      width=5, textvariable=self.EVENT)
        sb1.grid(row=3,column=3, sticky=W)
        sb1_label = Label(self, text="EVENT")
        sb1_label.grid(row=2,column=3, padx=0, sticky=W)
        
        
        
        self.LIMIT_L.set("19000")
        sb1 = Spinbox(self, from_=0, to=100000, 
                      width=5, textvariable=self.LIMIT_L)
        sb1.grid(row=5,column=2, sticky=W)
        sb1_label = Label(self, text="ROI Start ")
        sb1_label.grid(row=4,column=2, padx=0, sticky=W)
        
        self.LIMIT_H.set("22500")
        sb1 = Spinbox(self, from_=0, to=100000, 
                      width=5, textvariable=self.LIMIT_H)
        sb1.grid(row=5,column=3, sticky=W)
        sb1_label = Label(self, text="ROI End ")
        sb1_label.grid(row=4,column=3, padx=0, sticky=W)
        
        
        
           
        self.PULSE_R.set("20142")
        sb1 = Spinbox(self, from_=0, to=100000, 
                      width=8, textvariable=self.PULSE_R)
        sb1.grid(row=5,column=4, sticky=E)
        sb1_label = Label(self, text=" Pulse Rise")
        sb1_label.grid(row=4,column=4, padx=0, sticky=E)
        
        self.PULSE_L.set("1200")
        sb1 = Spinbox(self, from_=0, to=5000, 
                      width=8, textvariable=self.PULSE_L)
        sb1.grid(row=5,column=5, sticky=E)
        sb1_label = Label(self, text=" Pulse Length")
        sb1_label.grid(row=4,column=5, padx=0, sticky=E)
        
        
        
        sb1_label = Label(self, text="  ")
        sb1_label.grid(row=2,column=7, padx=0, sticky=W)        
        sb1_label = Label(self, text="  ")
        sb1_label.grid(row=6,column=7, padx=0, sticky=W)
        
        
        self.pulse_height.set("545.5")
        sb1 = Entry(self, width=8, textvariable=self.pulse_height)
        sb1.grid(row=7,column=3, sticky=E)
        sb1_label = Label(self, text=" Amplitude")
        sb1_label.grid(row=6,column=3, padx=0, sticky=E)
        
        self.amplitude_range.set("2")
        sb1 = Entry(self, width=8, textvariable=self.amplitude_range)
        sb1.grid(row=7,column=4, sticky=E)
        sb1_label = Label(self, text=" Loop Range")
        sb1_label.grid(row=6,column=4, padx=0, sticky=E)
        
        self.delta.set("0.1")
        sb1 = Entry(self, width=8, textvariable=self.delta)
        sb1.grid(row=7,column=5, sticky=E)
        sb1_label = Label(self, text=" Loop Delta")
        sb1_label.grid(row=6,column=5, padx=0, sticky=E)
        
        self.noise_sigma.set("4")
        sb1 = Entry(self, width=3, textvariable=self.noise_sigma)
        sb1.grid(row=5,column=6, sticky=E)
        sb1_label = Label(self, text=" Noise Threshold")
        sb1_label.grid(row=4,column=6, padx=0, sticky=E)
        
        sb_coeff_label = Label(self, text= "Coefficient ")
        sb_coeff_label.grid(row=0,column=6, padx=0, sticky=E)
        self.sb_coeff = Label(self)
        self.sb_coeff.grid(row=1,column=6, padx=0, sticky=E)


        
        # MAIN BUTTONS
        obtn = Button(self, text="GO!!", command=self.find_C)
        obtn.grid(row=14, column=4, sticky=E, pady=10)
        
        cbtn = Button(self, text="Quit", command=self.quit)
        cbtn.grid(row=14, column=5, sticky=E, pady=10)
        
        hbtn = Button(self, text="Help", command=self.help_f)
        hbtn.grid(row=14, column=0, sticky=W, pady=10)


    def help_f(self):
		top = Toplevel()
		top.title("HELP")
		msg = Message(top, width= 500,
             text="COEFF Calibration Procedure: \n \
             Input Start Point and Length of the pulse \n \
             Input an initial guess of the pulse amplitude \n \
             Use a ROI with at least 1000 samples of baseline \n \
             and 1000 samples after pulse end \n \
             Adjust loop range and step until a graph error \n \
             with a minimum is observed \n \
             Refine the search to increase precision")
		msg.pack()

		button = Button(top, text="Close", command=top.destroy)
		button.pack()

        
    
    def find_C(self):
        
        draw = False
        LIMIT_L       = self.LIMIT_L.get()      #19000
        LIMIT_H       = self.LIMIT_H.get()      #22500
        PULSE_R       = self.PULSE_R.get()      #20142
        PULSE_L       = self.PULSE_L.get()      #1200
        pulse_height  = self.pulse_height.get() #545
        hdf5_file 	 = self.hdf5_file.get()   #'2052.h5.z'
        PMT 		    = self.PMT.get()
        event         = self.EVENT.get()
        amplitude_range = self.amplitude_range.get() #2
        delta          = self.delta.get()       #0.1
        noise_sigma     = self.noise_sigma.get() #4
        
        
        coeff_aux = fc.find_coeff(LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
                             pulse_height,
                             hdf5_file, PMT, event,
                             amplitude_range, delta,
                             noise_sigma,
                             draw)
        
        plt.show() 
        self.sb_coeff.configure(text=str(coeff_aux))
																		
    
    
    def centerUI(self,w,h):
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()
        x  = (sw-w)/2
        y  = (sh-h)/2
        self.parent.geometry('%dx%d+%d+%d' % (w,h,x,y))
        
        

def main():

	root = Tk()
	root.resizable(width=False, height=False)
	app = gui(root)
	root.mainloop()


if __name__ == '__main__':
    main()        