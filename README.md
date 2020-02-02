# [Isochronous-Cyclotron](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/Isochronous_cyclotron.py)

```
Cyclotron and isochronous cyclotron particle accelerator sumilations. Implemented object-oriented 
methodology to resemblesubatomic particles and utilized the Runge-Kutta 4th order integration 
algorithm to solve electromagnetic differential equations and gather data of interest. Dynamic 
attribute updates allowed the model to accurately resemble relativistic principles and effects. 
```


```
A cyclotron, invented by Ernest O. Lawrence in 1934, is one of the earliest types of particle 
accelerators that utilizes electromagnetic fields to accelerate charged particles to extreme 
velocities [1]. A static magnetic field is applied inside two ”D-shaped” regions called dees,
which serves to keep the particle in a semi-circular path due to Lorentz Force.
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://www.codecogs.com/eqnedit.php?latex=\LARGE&space;F&space;=&space;qE&space;&plus;&space;qv&space;\times&space;B" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;F&space;=&space;qE&space;&plus;&space;qv&space;\times&space;B" title="\LARGE F = qE + qv \times B" /></a>

```
As the particle reaches a gap between the dees, it is accelerated by a rapidly varying electric 
field that increases both the particle’s velocity and radius. Each time the particle crosses a 
gap, the polarity of the electric field reverses to accelerate the particle in the correct direction. 
In order to properly time the polarity reversal, it is necessary for the electric field to be tuned 
to the cyclotron resonance. By equating the magnetic Lorentz force with the centripetal force and 
making the proper substitutions, one obtains the gyrofrequency
```
<img src="https://latex.codecogs.com/gif.latex?\LARGE&space;f=\frac{qB}{2\pi&space;m}" title="\LARGE f=\frac{qB}{2\pi m}" />

```
where q, B and m is the particle-charge, magnetic-field strength, and particle-mass, respectively. 
It is important to note that the gyrofrequency is independent of both the radius and velocity, 
and is thus contant over a static magnetic field. The resulting particle trajectory is an outward
spiral starting from the center of the cyclotron, as shown in figure 1.
```
![alt text](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/Particle_Trajec.png "Trajectory")

```
The cyclotron accelerator shown above utilizes a magnetic field with strength |B| = 1.5 T, pointing
along the y-direction and a square-wave electric field with strength |E⃗| = 5,000,000 N·C−1 in the
x-direction. The proton is injected at the center with initial velocity equal to 5% the speed of 
light and begins its spiral trajectory at a gyrofrequency equal to f = 22.87 MHz. The proton crosses 
the gap 46 times, that is to say, the iproton velocity (and hence energy) is increased a total number 
of 46 times before it is ejected. At the final jump, the particle is ejected with a final velocity of 
vf = 50, 196, 499 m/s — 16.7 % the speed of light — in a time of t = 1.09 μs, as shown in figure 2. 
The particle trajectory and velocity were computed using the Runge-Kutta 4th order numerical technique.
```

![alt text](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/velocity_vs_time.png "Velocity vs Time")

```
The gyrofrequency may also be obtained from experimental data from the particle, since velocity is 
proportional to the gyrofrequency and radius in the following manner
```
<img src="https://latex.codecogs.com/gif.latex?\LARGE&space;v=(2\pi&space;f)r" title="\LARGE v=(2\pi f)r" />
```
Therefore, plotting a graph of velocity as a function of radius gives a straight line with slope of 2πf. 
This allows a measurement of cyclotron frequency directly from real data. This plot is demonstrated in figure 3.
```

![alt text](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/veocity_vs_radius.png "Velocity vs. Rotation Radius")


# The Runge-Kutta 4th Order Integration Method
<details>
  <summary>Expand</summary>
  
```python
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
  ```
</details>

## References

* [1] U.S. Patent 1,948,384 Lawrence, Ernest O. Method and apparatus for the acceleration of ions, filed: January 26, 1932, granted: February 20, 1934
