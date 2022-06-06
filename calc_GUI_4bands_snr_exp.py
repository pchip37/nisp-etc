import tkinter 
from tkinter import *
from tkinter.ttk import Combobox
#from exptime_calc_NISP_v5_final_additions_imaging import *
from mag_exp_NISP_v5_final import *
from tkinter import messagebox
window = Tk()

#TITLE OF THE GUI WINDOW
window.title("NISP Exposure TIme Calculator (NISP ETC)")
window.geometry("1000x700+10+20")

#VARIOUS WIDGETS' PACKS ARE USED TO SHOW THE OBJECT IN THE GUI WINDOW



def welcome():
	wel = messagebox.showinfo("Hey!","Welcome to the NISP's Exposure time calculator!")

 	
	

#BUTTON WITH A CLICK ACTION
wel = Button(window, text = "Hello!",command = welcome)
wel.pack()
wel.place(x=450, y=50)


expos = Button(window, text = "Calculate",command = exposure_plot)
expos.pack()
expos.place(x=450, y=350)





def close_window():
	exit()
	
close = Button(window, text = "Close", command = close_window)
close.pack()
close.place(x=450, y=550)

window.mainloop()
