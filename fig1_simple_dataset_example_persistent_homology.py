import numpy as np
import matplotlib.pyplot as plt
import homcloud.interface as hc

data=([-0.5,0],[0,0],[0,1],[1,0],[1,1],[10,10],[10,15],[15,10],[15,15])
data=np.asarray(data)


r1=0.5
r2=np.sqrt(2)/2
r3=2.5
r4=5*np.sqrt(2)/2

plt.subplot(2,3,1)
plt.scatter(data[:,0],data[:,1],color='black')
plt.title('A) Example Data')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
 
plt.subplot(2,3,2)
plt.scatter(data[:,0],data[:,1],color='red')
c1=plt.Circle((0,0),radius=r1)
c2=plt.Circle((0,1),radius=r1)
c3=plt.Circle((1,0),radius=r1)
c4=plt.Circle((1,1),radius=r1)
c5=plt.Circle((10,10),radius=r1)
c6=plt.Circle((10,15),radius=r1)
c7=plt.Circle((15,10),radius=r1)
c8=plt.Circle((15,15),radius=r1)
c9=plt.Circle((-0.5,0),radius=r1)
plt.gca().add_artist(c1)
plt.gca().add_artist(c2)
plt.gca().add_artist(c3)
plt.gca().add_artist(c4)
plt.gca().add_artist(c5)
plt.gca().add_artist(c6)
plt.gca().add_artist(c7)
plt.gca().add_artist(c8)
plt.gca().add_artist(c9)
plt.grid(True)
plt.title('B) Radius '+str(r1))

plt.subplot(2,3,3)
plt.scatter(data[:,0],data[:,1],color='red')
c1=plt.Circle((0,0),radius=r2)
c2=plt.Circle((0,1),radius=r2)
c3=plt.Circle((1,0),radius=r2)
c4=plt.Circle((1,1),radius=r2)
c5=plt.Circle((10,10),radius=r2)
c6=plt.Circle((10,15),radius=r2)
c7=plt.Circle((15,10),radius=r2)
c8=plt.Circle((15,15),radius=r2)
c9=plt.Circle((-0.5,0),radius=r2)
plt.gca().add_artist(c1)
plt.gca().add_artist(c2)
plt.gca().add_artist(c3)
plt.gca().add_artist(c4)
plt.gca().add_artist(c5)
plt.gca().add_artist(c6)
plt.gca().add_artist(c7)
plt.gca().add_artist(c8)
plt.gca().add_artist(c9)
plt.grid(True)
plt.title('C) Radius '+str(round(r2,3)))

plt.subplot(2,3,4)
plt.scatter(data[:,0],data[:,1],color='red')
c1=plt.Circle((0,0),radius=r3)
c2=plt.Circle((0,1),radius=r3)
c3=plt.Circle((1,0),radius=r3)
c4=plt.Circle((1,1),radius=r3)
c5=plt.Circle((10,10),radius=r3)
c6=plt.Circle((10,15),radius=r3)
c7=plt.Circle((15,10),radius=r3)
c8=plt.Circle((15,15),radius=r3)
c9=plt.Circle((-0.5,0),radius=r3)
plt.gca().add_artist(c1)
plt.gca().add_artist(c2)
plt.gca().add_artist(c3)
plt.gca().add_artist(c4)
plt.gca().add_artist(c5)
plt.gca().add_artist(c6)
plt.gca().add_artist(c7)
plt.gca().add_artist(c8)
plt.gca().add_artist(c9)
plt.grid(True)
plt.title('D) Radius '+str(r3))

plt.subplot(2,3,5)
plt.scatter(data[:,0],data[:,1],color='red')
c1=plt.Circle((0,0),radius=r4)
c2=plt.Circle((0,1),radius=r4)
c3=plt.Circle((1,0),radius=r4)
c4=plt.Circle((1,1),radius=r4)
c5=plt.Circle((10,10),radius=r4)
c6=plt.Circle((10,15),radius=r4)
c7=plt.Circle((15,10),radius=r4)
c8=plt.Circle((15,15),radius=r4)
c9=plt.Circle((-0.5,0),radius=r4)
plt.gca().add_artist(c1)
plt.gca().add_artist(c2)
plt.gca().add_artist(c3)
plt.gca().add_artist(c4)
plt.gca().add_artist(c5)
plt.gca().add_artist(c6)
plt.gca().add_artist(c7)
plt.gca().add_artist(c8)
plt.gca().add_artist(c9)
plt.grid(True)
plt.title('E) Radius '+str(round(r4,3)))

plt.subplot(2,3,6)
output=hc.PDList.from_alpha_filtration(data,no_squared=True,save_boundary_map=True,save_phtrees=True,save_to="pointcloud.pdgm")
pd1=hc.PDList("pointcloud.pdgm").dth_diagram(1)
persist=np.asarray(pd1.deaths-pd1.births)
births=np.asarray(pd1.births)
deaths=np.asarray(pd1.deaths)
x=np.linspace(0,np.amax(births),100)
plt.scatter(pd1.births,pd1.deaths,color='blue',label='H1 features')
plt.scatter(x,x,color='green',label='y=x line')
plt.xlabel('Birth Radius')
plt.ylabel('Death Radius')
plt.legend()
plt.title('F) Dimension One Homology Features')

plt.tight_layout()
plt.show()

