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
![alt text](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/images/ForceEq.gif "Force")
```
As the particle reaches a gap between the dees, it is accelerated by a rapidly varying electric 
field that increases both the particle’s velocity and radius. Each time the particle crosses a 
gap, the polarity of the electric field reverses to accelerate the particle in the correct direction. 
In order to properly time the polarity reversal, it is necessary for the electric field to be tuned 
to the cyclotron resonance. By equating the magnetic Lorentz force with the centripetal force and 
making the proper substitutions, one obtains the gyrofrequency


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


Therefore, plotting a graph of velocity as a function of radius gives a straight line with slope of 2πf. 
This allows a measurement of cyclotron frequency directly from real data. This plot is demonstrated in figure 3.
```

![alt text](https://github.com/TonyAlarcon/Isochronous-Cyclotron/blob/master/veocity_vs_radius.png "Velocity vs. Rotation Radius")
