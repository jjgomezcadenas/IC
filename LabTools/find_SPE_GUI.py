# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 00:31:32 2016

@author: viherbos
"""


from PIL import Image, ImageTk
from Tkinter import Tk, BOTH, W, N, E, S
from Tkinter import Checkbutton, BooleanVar, IntVar, DoubleVar, Spinbox
from Tkinter import StringVar
from ttk import Frame, Button, Label, Entry
import find_SPE_h5 as SPE
import matplotlib.pyplot as plt

a=0

class Example(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
		self.initUI()
		self.centerUI(w=520,h=270)

	def initUI(self):

		self.parent.title("FIND SPE VALUE")
		self.pack(fill=BOTH, expand=True)

		self.columnconfigure(0, weight=1)
		#self.rowconfigure(0, weight=1)
		# weight attibute is used to make them growable

		self.graph_cb  = BooleanVar()
		self.bins      = IntVar()
		self.path      = StringVar()
		self.n_files   = IntVar()
		self.start_s   = IntVar()
		self.end_s     = IntVar()
		self.guess     = IntVar()


		search = Image.open("next_logo.jpg")
		search_temp = search.resize((160, 200), Image.ANTIALIAS)
		search_aux = ImageTk.PhotoImage(search_temp)
		label1 = Label(self, image=search_aux)
		label1.image = search_aux
		label1.grid(row=0, column=0,
				columnspan=10, rowspan=10, sticky=E+W+S+N)


		#Number of Files and Bins. Spin Box
		self.n_files.set("2000")
		sb1 = Spinbox(self, from_=1, to=10000,
				  width=6, textvariable=self.n_files)
		sb1.grid(row=1,column=4, sticky=W)
		sb1_label = Label(self, text="Files")
		sb1_label.grid(row=1,column=3, padx=5, sticky=E)

		self.bins.set("50")
		sb2 = Spinbox(self, from_=10, to=200,
				  width=6, textvariable=self.bins)
		sb2.grid(row=1,column=6, sticky=W)
		sb2_label = Label(self, text="Hist. Bins")
		sb2_label.grid(row=1,column=5, padx=5, sticky=E)

		# INTEGRATION LIMITS
		Integration_label = Label(self, text="INTEGRATION",
		                          font = "Verdana 12 bold")
		Integration_label.grid(row=3,column=4,
						padx=5,
						columnspan = 2)
		self.start_s.set("1732")
		sb3 = Spinbox(self, from_=1, to=10000,
				  width=6, textvariable=self.start_s)
		sb3.grid(row=4,column=4, sticky=W)
		sb3_label = Label(self, text="StartPoint")
		sb3_label.grid(row=4,column=3, padx=5, sticky=E)

		self.end_s.set("1752")
		sb4 = Spinbox(self, from_=1, to=10000,
				  width=6, textvariable=self.end_s)
		sb4.grid(row=4,column=6, sticky=W)
		sb4_label = Label(self, text="EndPoint")
		sb4_label.grid(row=4,column=5, padx=5, sticky=E)
		sb4_label = Label(self, text="")
		sb4_label.grid(row=4,column=7, padx=5, sticky=E)

		# FITTING PARAMETERS
		Integration_label = Label(self, text="FITTING",
		                          font = "Verdana 12 bold")
		Integration_label.grid(row=6,column=4,
						padx=5,
						columnspan = 2)
		self.guess.set("-20")
		sb5 = Spinbox(self, from_=-50, to=-1,
				  width=6, textvariable=self.guess)
		sb5.grid(row=7,column=4, sticky=W)
		sb5_label = Label(self, text="SPE guess")
		sb5_label.grid(row=7,column=5, padx=5, sticky=W)

		#Check buttons
		cb1 = Checkbutton(self, text="MultiGraph Output", variable=self.graph_cb					)
		cb1.select()
		cb1.grid(row=7,column=6, sticky=W)


		#Text Box
		#self.path.set("F:/DATOS_DAC/spe_1230/2046/pmt_0_trace_evt_")
		self.path.set("spe_1230_2046.h5.z")
		e1 = Entry(self, textvariable=self.path, width=45)
		e1.grid(row=10,column=3, sticky=W, columnspan=10, padx=10, pady=5)
		e1_label = Label(self, text="DataSet path (including name file)")
		e1_label.grid(row=9,column=3,sticky=W, columnspan=10, padx=10)


		# Main buttons
		obtn = Button(self, text="GO!!", command=self.SPE_f)
		obtn.grid(row=14, column=5, sticky=E, pady=5)

		cbtn = Button(self, text="Quit", command=self.quit)
		cbtn.grid(row=14, column=6, sticky=E, pady=5)

		hbtn = Button(self, text="Help")
		hbtn.grid(row=14, column=0, sticky=W, pady=5)

	def SPE_f(self):

		global a

		if (self.graph_cb.get()==True):
			b=a
			a=a+1
		else:
			b=0

		SPE.find_SPE_h5(base_path=self.path.get(),
				  n_files=self.n_files.get(),
				  start=self.start_s.get(),
				  end=self.end_s.get(),
				  bins=self.bins.get(),
				  guess=self.guess.get(),
				  n_figure=b)


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