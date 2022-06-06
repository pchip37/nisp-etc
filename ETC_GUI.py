################################################################################################################################################################################################################
######################### GUI for a Single-lined optical design with revised (based on manufacturer's responses on materials' availability) materials #############################################################
################################################################################################################################################################################################################
##################################################################################### Prachi Prajapati, March 10, 2022 ########################################################################################
################################################################################################################################################################################################################

import numpy as np
import numpy as geek
import scipy as scp
from scipy.constants import c,h,k
from scipy.integrate import quad
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from numpy import genfromtxt



#from exptime_calc_NISP_v5_final_additions_imaging import *
#from th_bg_temp_yjhkks_v1_filterwise_final import *



import tkinter 
from tkinter import ttk
from tkinter import *
from tkinter.ttk import Combobox
from tkinter import messagebox
window = Tk()

#TITLE OF THE GUI WINDOW
window.title("NISP Exposure TIme Calculator")
window.geometry("500x400+10+20")

#welcome
def welcome():
	wel = messagebox.showinfo("Hey!","Welcome to the NISP's Exposure time calculator!")    		

wel = Button(window, text = "Welcome!", bg="blue", fg="white", command = welcome)
wel.pack()
wel.place(x=225, y=5)

# Create text widget and specify size.
T = tkinter.Text(window, height = 5, width = 52)

#LABEL WIDGET
label_widget = tkinter.Label(window, text="User inputs", fg='green', font=("Helvetica", 16))
label_widget.pack() 
label_widget.place(x=213.5,y=40)

#Fact = """A man can be arrested in Italy for wearing a skirt in public."""

#BUTTON STYLES FOR THE GUI WINDOW
styleb = ttk.Style()

styleb.configure('W.TButton', font =
               ('calibri', 14, 'bold'),	#, 'underline'),
                foreground = 'blue')



#VARIOUS WIDGETS' PACKS ARE USED TO SHOW THE OBJECT IN THE GUI WINDOW

source_var = tkinter.StringVar()
temp_var = tkinter.DoubleVar()
snr_var = tkinter.DoubleVar()
mag_var = tkinter.DoubleVar()
filters_var = tkinter.StringVar()
za_var = tkinter.DoubleVar()

def submit_source():
	global source 
	
def submit_temp():
	global T 
	
def submit_snr():
	global SNR 
				
def submit_mag():
	global mag_obj_raw 
		
def submit_filters():
	global filters 
		
def submit_za():
	global ZA 




global n
n = 5.0		
#f-number based on camera. The telescope is f/8 though!
global D
D = 2.5	
#telescope primary diameter in m
global f
f = n*D		
#effective focal length in m
global D2
D2 = 1.002	
#telescope secondary diameter in m

global pixelsize
pixelsize = 18e-3		
#in mm

global platescale
#platescale = 206265./(f*10**3)*pixelsize	
platescale = 0.29 #FOV(x)/PIXELS(x) = 600./2048.
#in arcsec/pixel





	
def submit():
	source = source_var.get()
	T = float(temp_var.get())
	SNR = float(snr_var.get())
	mag_obj_raw = float(mag_var.get())
	filters = filters_var.get()
	ZA = float(za_var.get())
	
	print("Name of your source is : "+source)
	source_var.set("")
	print("Ambient tempearture is : ",str(T))
	temp_var.set("")
	print("S/N ratio is : "+str(SNR))
	snr_var.set("")	
	print("Magnitude of the object of interest is : "+str(mag_obj_raw))
	mag_var.set("")
	print("Chosen filter is : "+filters)
	filters_var.set("")
	print("Current zenith angle of the object is : "+str(ZA))
	za_var.set("")




	#tf is the lens %T etc.
	#tf = MR_1*MR_2*coll_lenses*cam_lenses*foldmirr
	#commenting this because the view-factor takes care of all these parameters
	global tf 
	tf = 1.0
	global vf #view-factor with a pupil-mask        
	vf = 4.6*(10**-7)
	#vf = 3.27*(10**-6) #view-factor without a pupil-mask



	if filters in ['Y','y']:
		T = [T]
#		global lambdaa  #central wavelength in microns
		lambdaa = 1.02
#		global dellambda #bandwidth in microns
		dellambda = 0.25
#		global QE  #quantum efficiency of detector (counts/photos)
		QE = 0.86 
#		global scale_factor  #in W/sq.m./nm
		scale_factor = 5.71e-12 
#		global s  #sky brightness in mag/sq.arcsec
		s = 17.0
#		global extinction  #extinction co-efficient with units in magnitude/airmass (For Abu, Rajasthan) - avg. value, J -band (Gupta et al., 2004) - other bands are extended from it and some taken from KPNO data (since the J-value was closeby)
		extinction = 0.2
#		global fwhm  #in pixels (seeing limited - taken the worst case 1.04" seeing)
		fwhm = 1.04*(platescale**-1)
#		global d 
		d = 0.001
#		global l	#wavelength in microns
		l = np.arange(0.952,1.091,d)
		
#		global t 
		t = np.genfromtxt('Yband.csv', delimiter=",", skip_header=1, names="l,t", usecols=("t"))
		
		def bb(l,T):
			bb = ((2.*h*c**2)/((l*10**-6)**5))*((1.)/(np.exp((h*c)/((l*10**-6)*k*T))-(1.)))*(10**-6)
			return bb  #Planck's function in W/sq.m./micron/sr		
		def convol(l,T,t):
			convol = bb(l,T)*t["t"]/100.
		#	convol = np.multiply((bb(lambdaa,T)),th)
			return convol		
#		global S
		S = []
#		global ss
		ss = np.zeros((len(T),len(l)))		
		for z in range(len(T)):
			for p in range(len(l)):
				ss[z][p] = ((((convol(l,T[z],t)[p])*d)))
		#	print(sum(syy[z]))
			S.append(float(sum(ss[z])))	 
		S = np.array(S)
		#print(len(s))		
		#effective telescope structural area contributing to thermal radiation is considered to be 5 sq.m.
#		global s_px #counts/sec per pixel
		s_px = S*QE*((tf*np.pi*5.*vf)/((2048.)**2))/((h*c/(lambdaa*10**-6)))
#		global th  #thermal counts/sec/pixel
		th = s_px[0]
	
	elif filters in ['J','j']:
		T = [T]	
		lambdaa = 1.25
		dellambda = 0.16
		QE = 0.87
		scale_factor = 2.9e-12
		s = 15.13
		extinction = 0.15
		fwhm = 0.999*(platescale**-1)
		d = 0.001
		l = np.arange(1.147,1.369,d)
		t = np.genfromtxt('Jband.csv', delimiter=",", skip_header=1, names="l,t", usecols=("t"))
		def bb(l,T):
			bb = ((2.*h*c**2)/((l*10**-6)**5))*((1.)/(np.exp((h*c)/((l*10**-6)*k*T))-(1.)))*(10**-6)
			return bb  #Planck's function in W/sq.m./micron/sr		
		def convol(l,T,t):
			convol = bb(l,T)*t["t"]/100.
			return convol		
		S = []
		ss = np.zeros((len(T),len(l)))		
		for z in range(len(T)):
			for p in range(len(l)):
				ss[z][p] = ((((convol(l,T[z],t)[p])*d)))
			S.append(float(sum(ss[z])))	 
		S = np.array(S)
		s_px = S*QE*((tf*np.pi*5.*vf)/((2048.)**2))/((h*c/(lambdaa*10**-6)))
		th = s_px[0]
	
	elif filters in ['H','h']:
		T = [T]
		lambdaa = 1.64
		dellambda = 0.29
		QE = 0.84
		scale_factor = 1.1e-12 
		s = 13.66
		extinction = 0.09
		fwhm = 0.946*(platescale**-1)
		d = 0.0001
		l = np.arange(1.4622,1.8261,d)
		t = np.genfromtxt('Hband.csv', delimiter=",", skip_header=1, names="l,t", usecols=("t"))
		def bb(l,T):
			bb = ((2.*h*c**2)/((l*10**-6)**5))*((1.)/(np.exp((h*c)/((l*10**-6)*k*T))-(1.)))*(10**-6)
			return bb  #Planck's function in W/sq.m./micron/sr		
		def convol(l,T,t):
			convol = bb(l,T)*t["t"]/100.
			return convol		
		S = []
		ss = np.zeros((len(T),len(l)))		
		for z in range(len(T)):
			for p in range(len(l)):
				ss[z][p] = ((((convol(l,T[z],t)[p])*d)))
			S.append(float(sum(ss[z])))	 
		S = np.array(S)			
		s_px = S*QE*((tf*np.pi*5.*vf)/((2048.)**2))/((h*c/(lambdaa*10**-6)))
		th = s_px[0] 
		
	elif filters in ['KS','ks','Ks']:
		T = [T]
		lambdaa = 2.15
		dellambda = 0.32
		QE = 0.87
		scale_factor = 4.18e-13 
		s = 12.28
		extinction = 0.074
		fwhm = 0.896*(platescale**-1)
		d = 0.001
		l = np.arange(1.944,2.374,d)
		t = np.genfromtxt('Ksband.csv', delimiter=",", skip_header=1, names="l,t", usecols=("t"))
		def bb(l,T):
			bb = ((2.*h*c**2)/((l*10**-6)**5))*((1.)/(np.exp((h*c)/((l*10**-6)*k*T))-(1.)))*(10**-6)
			return bb  #Planck's function in W/sq.m./micron/sr		
		def convol(l,T,t):
			convol = bb(l,T)*t["t"]/100.
			return convol		
		S = []
		ss = np.zeros((len(T),len(l)))		
		for z in range(len(T)):
			for p in range(len(l)):
				ss[z][p] = ((((convol(l,T[z],t)[p])*d)))
			S.append(float(sum(ss[z])))	 
		S = np.array(S)
		s_px = S*QE*((tf*np.pi*5.*vf)/((2048.)**2))/((h*c/(lambdaa*10**-6)))
		th = s_px[0]
		
	elif filters in ['K','k']:
		T = [T]
		lambdaa = 2.20 
		dellambda = 0.34
		QE = 0.87
		scale_factor = 4.18e-13 
		s = 12.28
		extinction = 0.074
		fwhm = 0.896*(platescale**-1)
		d = 0.0005
		l = np.arange(1.9875,2.4005,d)
		t = np.genfromtxt('Kband.csv', delimiter=",", skip_header=1, names="l,t", usecols=("t"))
		def bb(l,T):
			bb = ((2.*h*c**2)/((l*10**-6)**5))*((1.)/(np.exp((h*c)/((l*10**-6)*k*T))-(1.)))*(10**-6)
			return bb  #Planck's function in W/sq.m./micron/sr		
		def convol(l,T,t):
			convol = bb(l,T)*t["t"]/100.
			return convol	
		S = []
		ss = np.zeros((len(T),len(l)))		
		for z in range(len(T)):
			for p in range(len(l)):
				ss[z][p] = ((((convol(l,T[z],t)[p])*d)))
			S.append(float(sum(ss[z])))	 
		S = np.array(S)								
		s_px = S*QE*((tf*np.pi*5.*vf)/((2048.)**2))/((h*c/(lambdaa*10**-6)))
		th = s_px[0]
		
	else:
		print("Choose the correct filter.")	


	p = 1.4*(fwhm**2)	
	#pixel area covered by a point-source : units in (pixels)**2

	global rn 
	rn = 7./np.sqrt(p)
	#read-noise defined as [sqrt(p*R**2)] where, R is read-noise (in electrons/pixel) @Fowler-8 sampling
	global binning
	binning = 1.00 
	#binning factor

	
	#######################################################################
	
	#transmission (0.24) through collimator lenses - total eight lenses
	global zns_coll 
	zns_coll = 0.75**4
	global caf2_coll 
	caf2_coll = 0.97**1
	global fused_silica_coll
	fused_silica_coll = 0.98**1
	global baf2_coll
	baf2_coll = 0.9**2


	#transmission (0.94) through three filters 
	global fused_silica_filter 
	fused_silica_filter = 0.98**3


	#transmission (0.33) through camera lenses - total six lenses
	global baf2_cam 
	baf2_cam = 0.9**2
	global zns_cam
	zns_cam = 0.75**3
	global fused_silica_cam 
	fused_silica_cam = 0.98**1


	#primary mirror reflectivity
	global MR_1
	MR_1 = 0.95

	#secondary mirror reflectivity
	global MR_2 
	MR_2 = 0.95 		

	#secondary mirror unvignetting factor - dep. on how much light it allows
	global unvig_factor
	unvig_factor = 0.75


	global coll_lenses 
	coll_lenses = 	zns_coll*caf2_coll*fused_silica_coll*baf2_coll

	global Filters
	Filters = fused_silica_filter		

	global cam_lenses 
	cam_lenses = baf2_cam*zns_cam*fused_silica_cam


	#efficiency of telescope -----> Total around 5.2%
	global tel_efficiency 
	tel_efficiency = MR_1*MR_2*unvig_factor*coll_lenses*cam_lenses*Filters


	global spiders_area 	#sq.m.
	spiders_area = 4.*(15*10**-3)*(0.749)	


	#telescope's effective light collecting area in sq.m. (considered the blockage by the secondary mirror)
	global tel_area 
	tel_area = (0.25*np.pi*(D**2-D2**2)) - (spiders_area)

	#airmass at given zenith-angle
	global X 	#*(1-(0.0012*((np.cos(ZA*np.pi/180.))**-2-1)))
	X = ((np.cos(ZA*np.pi/180.))**(-1))

	global m_obj_corr 
	m_obj_corr = mag_obj_raw+1.086*extinction*X

	global dark 	#dark current in photons/sec/pixel @77K
	dark = 0.01/QE
	
	global det_gain  #detector gain in photons/ADU		
	det_gain = 3.0/QE





	def F0(lambdaa,dellambda,scale_factor):
		F0 = (scale_factor*dellambda*10**3)*((h*c/(lambdaa*10**-6))**-1)
		return F0
		#flux density of A0 Vega like point-object in photons/s/sq.m. = scale_band*lambdaa_bandwidth in (microns)	
	
	def F_obj(lambdaa,dellambda,scale_factor,m_obj_corr):	
		F_obj = F0(lambdaa,dellambda,scale_factor)*(10**((-1.00*(m_obj_corr-0.00))/2.5))
		#flux density from object in photons/s/sq.m.
		return F_obj

	def N_obj(lambdaa,dellambda,scale_factor,m_obj_corr,tel_efficiency,tel_area,QE):	
		N_obj = F_obj(lambdaa,dellambda,scale_factor,m_obj_corr)*tel_efficiency*tel_area*QE
		return N_obj
		#object counts in counts/sec converted from counts/sec
	
	def SF_obj(lambdaa,dellambda,scale_factor,platescale,binning):
		SF_obj = F0(lambdaa,dellambda,scale_factor)*(10**((-1.00*(s-0.00))/2.5))*(platescale**2)*(binning**2)
		return SF_obj
		
	def S_obj(lambdaa,dellambda,scale_factor,platescale,binning,tel_efficiency,tel_area,QE):
		S_obj = SF_obj(lambdaa,dellambda,scale_factor,platescale,binning)*tel_efficiency*tel_area*QE
		return S_obj		
	
	
	def int_time(m_obj_corr,SNR,lambdaa,dellambda,scale_factor,platescale,binning,tel_efficiency,tel_area,QE,dark,th,p,rn,det_gain):
		A = (N_obj(lambdaa,dellambda,scale_factor,m_obj_corr,tel_efficiency,tel_area,QE))**2
		B = (-1*(SNR)**2)*(N_obj(lambdaa,dellambda,scale_factor,m_obj_corr,tel_efficiency,tel_area,QE)+	(S_obj(lambdaa,dellambda,scale_factor,platescale,binning,tel_efficiency,tel_area,QE)+th+dark)*p)	
		C = (-1*(SNR)**2)*(((p)*((rn)**2))+(p*(det_gain/2.)**2))
		t = np.abs(((-B)+(np.sqrt(((B)**2)-(4*A*C))))/(2*A))
		#required exposure time (in seconds)
		return t
	
	def printing():
		inttime = 	int_time(m_obj_corr,SNR,lambdaa,dellambda,scale_factor,platescale,binning,tel_efficiency,tel_area,QE,dark,th,p,rn,det_gain)
		inttime = "{:.2f}".format(inttime)
		return inttime
			
	
	global x 
	x = printing()
#	print("{:.2f}".format(x),'sec')	
#	tkinter.messagebox.showinfo("The required exposure time for "+source+" is ")
	tkinter.messagebox.showinfo("Result", "Your chosen S/N ratio is "+str(SNR)+". \n\nExpected ambient temperature is "+str(int(T[0]))+". \n\nIn "+filters+" filter, the required exposure time for "+source+" (with magnitude "+str(mag_obj_raw)+") at zenith angle of "+str((ZA))+" is "+str(x)+" seconds.")
	






			
	
source_label = ttk.Label(window, text = 'Source', font=('calibre',10, 'bold'))
source_entry = ttk.Entry(window,textvariable = source_var, font=('calibre',10,'normal'))

temp_label = ttk.Label(window, text = 'Ambient temperature', font=('calibre',10, 'bold'))
temp_entry = ttk.Entry(window,textvariable = temp_var, font=('calibre',10,'normal'))

snr_label = ttk.Label(window, text = 'S/N ratio', font=('calibre',10, 'bold'))
snr_entry = ttk.Entry(window,textvariable = snr_var, font=('calibre',10,'normal'))

mag_label = ttk.Label(window, text = 'Object magnitude', font=('calibre',10, 'bold'))
mag_entry = ttk.Entry(window,textvariable = mag_var, font=('calibre',10,'normal'))

filters_label = ttk.Label(window, text = 'Filter', font=('calibre',10, 'bold'))
filters_entry = ttk.Entry(window,textvariable = filters_var, font=('calibre',10,'normal'))

za_label = ttk.Label(window, text = 'Zenith angle', font=('calibre',10, 'bold'))
za_entry = ttk.Entry(window,textvariable = za_var, font=('calibre',10,'normal'))


sub_btn = Button(window,text = 'Submit', fg="yellow", bg="green", command = submit)


source_label.pack()
source_label.place(x=100, y=75)

source_entry.pack()
source_entry.place(x=300, y=75)


temp_label.pack()
temp_label.place(x=100, y=100)

temp_entry.pack()
temp_entry.place(x=300, y=100)


snr_label.pack()
snr_label.place(x=100, y=125)

snr_entry.pack()
snr_entry.place(x=300, y=125)


mag_label.pack()
mag_label.place(x=100, y=150)

mag_entry.pack()
mag_entry.place(x=300, y=150)


filters_label.pack()
filters_label.place(x=100, y=175)

filters_entry.pack()
filters_entry.place(x=300, y=175)


za_label.pack()
za_label.place(x=100, y=200)

za_entry.pack()
za_entry.place(x=300, y=200)


sub_btn.pack()
sub_btn.place(x=220, y=250)




#close
def close_window():
	exit()
	
close = Button(window, text = "Close", bg="blue", fg="white", command = close_window)
close.pack()
close.place(x=225, y=300)




#FOR A COMBO BOX WITH MULTIPLE OPTIONS - ONE/MULTIPLE TO BE SELECTED AMONG THESE OPTIONS
#var = StringVar()
#var.set("one")

'''
filter_label = ttk.Label(window, text = 'Filter', font=('calibre',10, 'bold'))
filter_var = ("Y", "J", "H", "Ks/K")
filters = Combobox(window, values=filter_var)
'''

#filter_entry = ttk.Entry(window,textvariable = filter_var, font=('calibre',10,'normal'))

'''
filter_label.place(x=55,y=300)
filters.place(x=100,y=300)	
'''	

window.mainloop()	
