import numpy as np
import math as math 
from cmath import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal

c = 3.0E08
b = [0.0, 0.0, -.3]
bmag = np.linalg.norm(b)

class Particle: 
	def __init__(self, pos, vel, mass, charge): 
	
		
		self.pos = pos
		self.vel = vel
		self.mass = mass
		self.charge = charge
	
	def gamma_factor(self):
		vel_mag = np.linalg.norm(np.array(self.vel))
		gamma = 1.0/ (sqrt(1 - (vel_mag/c)**2))
		return gamma
		
	def radius(self,b):
		b =  np.linalg.norm(b)
		speed = np.linalg.norm(self.vel)
		radius = speed * self.mass / (self.charge * b)
		return radius
	
	@property
	def period(self):
		field = np.linalg.norm(b)* self.gamma_factor() 
		period = 2.0 * math.pi * self.mass /( field * self.charge )
		return period
		
	@period.setter
	def period(self,mass,charge):
		self.mass = mass
		self.charge = charge
		self.vel
		
	



# choose a proton as the particle, describe in 3 Dimension 
proton = Particle([0.00, 0.0, 0.0], [0.05*c, 0.0, 0.0], 1.67E-27, +1.60E-19)

#print proton.period
gap_size = .5 * 0.104375 #proton.radius(b) = 0.104375 for B = -1.5 T

print proton.radius(b)
print proton.period


###############################################################################
#  
#   i determines which time-step are we on. For example if RK4 is in iteration of 
#	i = 200, then t = 200*h = half_a_period
#
###############################################################################


def E_field(particle,i): 
	V_0 = 1.5625e27 #Value that gives E = 2.5E8 when cos() = 1
	T = 4.37204977625e-08
	n = 400
	h = 4.37204977625e-08 / n
	t = np.zeros([60033])
	t[i] = i * h
	ang_frequency = ( 2 * math.pi)/ particle.period
	V = V_0 * np.cos(ang_frequency * t[i] ) 
	
	E = [V * 1.6e-19, 0 , 0] #Volts/meter
	
	return E
	

###################################################################################
#	Given position and velocity, the function returns acceleration as an array in 3D
#	The magnetic field is applied where x > 0 or  x < - gap_size
#	The electric field is applied where x < 0 and x > - gap_size
###################################################################################

jumps =0
jumps_max = 1000
count = 0

def a( q_over_m, r, v , i):
	global count, jumps                # if you use the global statement, the variable will become available 
	if jumps >= jumps_max:			# "outside" the scope of the function, effectively becoming a global variable
		# No acceleration
		a = 0.0
	elif r[0] >= 0 or r[0] <= -gap_size:
		a = np.cross(v, b) 
		a = a * q_over_m
		if count:
			jumps += 1
			velo[jumps] = np.linalg.norm(v)
			count = 0
	else: 	
		a = np.array(E_field(proton,i))
		a = a * q_over_m
		if r[1] > 0:
			a = a
		count += 1
	return a


###############################################################################
#  
#    Solving Equations of Motions Using RK4
#
###############################################################################



# initializes array that will hold the velocity of particle 
# each time it gains velocity (crosses gap)
velo = np.zeros((jumps_max + 1,))
velo[0] = np.linalg.norm(proton.vel)


def rk4(particle, iterations, desired_value):

	RK4_pos = []
	RK4_vel = []
	
	
	n = 400
	h = 4.37204977625e-08 / n #initial period
	q_over_m = 1.6e-19 / 1.67E-27
	
	p0 = np.array(particle.pos)
	v0 = np.array(particle.vel) 
	for i in range(iterations):
		#vel_mag = np.linalg.norm(np.array(particle.vel))
		#comment the following block to neglect relativistic effects
		#print 'gamma is', gamma
		particle.mass =  1.67E-27 #this is important otherwise gives wrong mass update
		particle.mass = particle.gamma_factor()*particle.mass
		#print 'new period is', particle.period
		#print 'new mass is',particle.mass
		#print 'particle velocity is', vel_mag/c *100, 'percent the speed of light' 

		p1 = p0 
		v1 = np.array(particle.vel)
		a1 = h * a( q_over_m, p1, np.array(particle.vel) , i )
		v1 = h * v1
		
		p2 = p0 + (v1 * 0.5)
		v2 = np.array(particle.vel) + (a1 * 0.5)
		a2 = h * a( q_over_m, p2, v2 , i)
		v2 = h * v2
		
		p3 = p0 + (v2 * 0.5)
		v3 = np.array(particle.vel) + (a2 * 0.5)
		a3 = h * a( q_over_m, p3, v3 , i)
		v3 = h * v3
		
		p4 = p0 + v3
		v4 = np.array(particle.vel) + a3
		a4 = h * a( q_over_m, p4, v4, i )
		v4 = h * v4
		
		particle.vel = np.array(particle.vel) + (a1 + 2.0 * (a2 + a3) + a4) / 6.0
		p0 = p0 + (v1 + 2.0 * (v2 + v3) + v4) / 6.0
	
		i += 1
		
		if desired_value == 'velocity': 
			vel_mag = np.linalg.norm(np.array(particle.vel))
			RK4_vel.append(vel_mag)
			value = RK4_vel
		elif desired_value == 'position': 
			RK4_pos.append(p0)
			value = RK4_pos
	
		
	return value
	

#######################################################################################
#  
#    Since this program graphs many figures, it is more convenient and efficient
#		to create a plot function. 
#		This function can accomodate 2D and 3D plots, note the first argument
#		must be a 3 element list, whose entries are labels for x,y,z coordinates, respectively.
#
#########################################################################################

def plot(coord_labels, title, dimension, x, y, z = 0):
	
	fig = plt.figure()
	
	if dimension == '2d':
		ax = fig.add_subplot(1,1,1)
		ax.plot(x, y )
	elif dimension == '3d':
		ax = fig.add_subplot(111,projection = dimension)
		ax.plot3D(x, z, y)
		ax.set_zlabel(coord_labels[2])
		
	ax.set_xlabel(coord_labels[0])
	ax.set_ylabel(coord_labels[1])
	ax.set_title(title)
	
		


###############################################################################
#  
#    Plotting Trajectory of Particle
#
###############################################################################


def position_plot(particle):

	x = []
	y = []
	z = []
	
	x.append(particle.pos[0])
 	y.append(particle.pos[1])
 	z.append(particle.pos[2])
	
	results = rk4(proton , 40032 ,'position')	
	
	for i in results:
		x.append(i[0])
		y.append(i[1])
		z.append(i[2])
		
	
	coordlabel = ["X-axis (m)","Y-axis (m)","Z-axis (m)"]
	plot(coordlabel, "Particle Trajectory in Cyclotron Accelerator", '3d', x,y,z)

	

###############################################################################
#  
#    velocity	vs. time plot
#
###############################################################################


def velocity_plot(particle):
	n = 400
	h = 4.37204977625e-08 / n
	
	v = rk4(proton , 60032 ,'velocity')	
	t = np.linspace(0, len(v) * h, len(v))
	speed_of_light_percentage = np.array(v)/c
	
	
	coordlabel = ["Time (s)","Velocity (m/s)","Z-axis (m)"]
	plot(coordlabel, "Beam Velocity as a Function of Time", '2d', t,speed_of_light_percentage)
	

	#prints final velocity
	final_velo = velo[len(velo) - 1]
 	print "final velocity is" ,final_velo/c

	#prints final time
	final_t = t[len(velo) - 1]
	print "final time is" ,final_t

###############################################################################
#  
#    velocity vs radius
#
###############################################################################


def vel_rad_plot(particle):


	radius = np.zeros((len(velo),))

	for item in range(1, len(velo)+1):
		radius[item-1] = velo[item-1] * particle.mass / (particle.charge * bmag)
	
	coordlabel = ["Radius (m)","Velocity (m/s)","Z-axis (m)"]
	plot(coordlabel, "Beam Velocity at Various Extraction Radius", '2d', radius,velo)
	

# 
# velocity_plot(proton)
# jumps = 0 
#position_plot(proton)
print rk4(proton , 40032 ,'position')[10]
# # vel_rad_plot(proton)
# plt.show()



