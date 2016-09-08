# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 01:33:41 2016

@author: viherbos
"""


from PIL import Image, ImageTk
from Tkinter import Tk, BOTH, W, N, E, S
from Tkinter import Checkbutton, BooleanVar, IntVar, DoubleVar, Spinbox
from Tkinter import StringVar, Toplevel, Message
from ttk import Frame, Button, Label, Entry
import DBLR as DB
from panel_to_hdf5 import read_panel_hdf5
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tkMessageBox

a=0

class Example(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
		self.initUI()
		self.centerUI(w=560,h=235)

	def initUI(self):

		self.parent.title("DBLR for Dummies")
		self.pack(fill=BOTH, expand=True)

		#self.columnconfigure(0, weight=1)
		#self.rowconfigure(0, weight=1)
		# weight attibute is used to make them growable

		self.meas       = IntVar()
		self.point      = IntVar()
		self.base_path  = StringVar()
		self.coef       = DoubleVar()
		self.noise      = DoubleVar()
		self.n_sigma    = DoubleVar()
		self.thr1       = DoubleVar()
		self.thr2       = DoubleVar()
		self.thr3       = DoubleVar()
		self.graph_sw   = BooleanVar()

		#Image
		factor=0.65
		search = Image.open("NEXT_official_logo.jpg")
		width_org, height_org = search.size
		search_temp = search.resize((int(width_org*factor),
							int(height_org*factor)), Image.ANTIALIAS)
		search_aux = ImageTk.PhotoImage(search_temp)
		label1 = Label(self, image=search_aux)
		label1.image = search_aux
		label1.grid(row=0, column=4,
				columnspan=2, rowspan=3, padx=5)


		#self.base_path.set("F:/DATOS_DAC/2052/pmt_0_trace_evt_")
		self.base_path.set("Argon.h5.z")
		e1 = Entry(self, textvariable=self.base_path, width=40)
		e1.grid(row=0,column=1, sticky=W, columnspan=3, pady=5, padx=5)
		e1_label = Label(self, text="Path & Name")
		e1_label.grid(row=0,column=0, columnspan=1, sticky=E, pady=5, padx=5)

		self.point.set("0")
		sb1 = Spinbox(self, from_=0, to=100,
				  width=4, textvariable=self.point)
		sb1.grid(row=2,column=1, sticky=W, pady=5, padx=5)
		sb1_label = Label(self, text="PMT")
		sb1_label.grid(row=2,column=0, padx=5, sticky=E)

		self.meas.set("0")
		sb1 = Spinbox(self, from_=0, to=100,
				  width=4, textvariable=self.meas)
		sb1.grid(row=2,column=3, sticky=W, pady=5, padx=5)
		sb1_label = Label(self, text="Event")
		sb1_label.grid(row=2,column=2, padx=5, sticky=E)

		#Check buttons
#		cb1 = Checkbutton(self, text="New Graph", variable=self.graph_sw)
#		cb1.select()
#		cb1.grid(row=2,column=2, sticky=W)


		#PARAMETERS
		Integration_label = Label(self, text="PARAMETERS",
		                          font = "Verdana 12 bold")
		Integration_label.grid(row=3,column=1,
						padx=5,
						columnspan = 2, pady=10, sticky=E)

		self.coef.set("1.65E-3")
		e2 = Entry(self, width=12, textvariable=self.coef)
		e2.grid(row=4,column=1, sticky=W, pady=5, padx=5)
		e2_label = Label(self, text="DBLR Coef")
		e2_label.grid(row=4,column=0, sticky=E, pady=5, padx=5)

		self.noise.set("0.75")
		e3 = Entry(self, width=12, textvariable=self.noise)
		e3.grid(row=4,column=3, sticky=W, pady=5, padx=5)
		e3_label = Label(self, text="Noise (LSB)")
		e3_label.grid(row=4,column=2, sticky=E, pady=5, padx=5)

		self.n_sigma.set("4")
		e4 = Entry(self, width=12, textvariable=self.n_sigma)
		e4.grid(row=4,column=5, sticky=W, pady=5, padx=5)
		e4_label = Label(self, text="Noise Threshold")
		e4_label.grid(row=4,column=4, sticky=E, pady=5, padx=5)

		self.thr1.set("0")
		e5 = Entry(self, width=12, textvariable=self.thr1)
		e5.grid(row=5,column=1, sticky=W, pady=5, padx=5)
		e5_label = Label(self, text="Threshold 1")
		e5_label.grid(row=5,column=0, sticky=E, pady=5, padx=5)

		self.thr2.set("0")
		e6 = Entry(self, width=12, textvariable=self.thr2)
		e6.grid(row=5,column=3, sticky=W, pady=5, padx=5)
		e6_label = Label(self, text="Threshold 2")
		e6_label.grid(row=5,column=2, sticky=E, pady=5, padx=5)

		self.thr3.set("0")
		e7 = Entry(self, width=12, textvariable=self.thr3)
		e7.grid(row=5,column=5, sticky=W, pady=5, padx=5)
		e7_label = Label(self, text="Threshold 3")
		e7_label.grid(row=5,column=4, sticky=E, pady=5, padx=5)





		# Main buttons
		obtn = Button(self, text="GO!!", command=self.DBLR_f)
		obtn.grid(row=6, column=4, sticky=E, pady=10)

		cbtn = Button(self, text="Quit", command=self.quit)
		cbtn.grid(row=6, column=5, sticky=E, pady=10)

		hbtn = Button(self, text="Help", command=self.help_f)
		hbtn.grid(row=6, column=0, sticky=W, pady=10)


	def help_f(self):
		top = Toplevel()
		top.title("HELP")
		msg = Message(top, width= 500,
                 text="Noise Threshold (NT) - Noise Based Main Threshold \
(Sigmas)\n Threshold 1 - Main Threshold (LSBs) \n Threshold 2 - End of Pulse Error \n \
Threshold 3 - End of Tail Error \n When Thr1 = Thr3 = 0 their values are defined as: \n \
(THR1 = NT (LSBs) / THR3 = NT*NOISE_ADC / 5)")
		msg.pack()

		button = Button(top, text="Close", command=top.destroy)
		button.pack()


	def DBLR_f(self):

		global a

		#if (self.graph_sw.get()==True):
		b=a
		a=a+1
		#else:
		#	b=0

		aux=self.base_path.get()
#		path=''.join([aux,str(self.meas.get()),'.txt'])
#		g=pd.read_csv(path)
		f=read_panel_hdf5(aux, self.point.get(), self.meas.get())
		f=4096-f

		recons,energia = DB.BLR( signal_daq=f.flatten().astype(float),
						  coef=self.coef.get(),
						  n_sigma=self.n_sigma.get(),
						  NOISE_ADC=self.noise.get(),
						  thr1=self.thr1.get(),
						  thr2=self.thr2.get(),
						  thr3=self.thr3.get(),
						  plot=False)
		plt.figure(a)
		plt.plot(recons)
		plt.figtext(0.2,0.75, ('ENERGY = %0.2f ' % (energia)))
		plt.show()


	def centerUI(self,w,h):
		sw = self.parent.winfo_screenwidth()
		sh = self.parent.winfo_screenheight()
		x  = (sw-w)/2
		y  = (sh-h)/2
		self.parent.geometry('%dx%d+%d+%d' % (w,h,x,y))

def main():

	root = Tk()
	root.resizable(width=False, height=False)
	app = Example(root)
	root.mainloop()


if __name__ == '__main__':
    main()