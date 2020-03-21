#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 22:03:48 2020

@author: bikramp
"""
import time as it

start = it.time()

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import os
#----------given file ----is -stud15----

zr=float(input('enter required redshift(0<z<=5)'))
#..........Beheroozi fuction for assigning stellar mass at redshift 5 and newly created hallos...........
class Behroozi_2013a:
            def __init__(self):
                self.eps0=-1.777;
                self.epsa=-0.006;
                self.epsz=-0.000;
                self.epsa2=-0.119;
        
                self.xM10=11.514;
                self.xM1a=-1.793;
                self.xM1z=-0.251;
        
                self.alp0=1.412;
                self.alpa=-0.731;
        
                self.delta0=3.508;
                self.deltaa=2.608;
                self.deltaz=-0.043;
        
                self.gamma0=0.316;
                self.gammaa=1.319;
                self.gammaz=0.279;
        
        
            def evolved_factors(self,zz):
                a=1./(1.+zz);
                nu=np.exp(-4*a*a);
                alp=self.alp0+self.alpa*(a-1.0)*nu;
                delta=self.delta0+(self.deltaa*(a-1.)+self.deltaz*zz)*nu;
                gamma=self.gamma0+(self.gammaa*(a-1.)+self.gammaz*zz)*nu;
                eps=self.eps0+(self.epsa*(a-1.0)+self.epsz*zz)*nu+self.epsa2*(a-1.)
                xM1=self.xM10+(self.xM1a*(a-1.0)+self.xM1z*zz)*nu;
                return nu, alp, delta, gamma, eps, xM1
        
            def fbeh12(self, x, nu, alp, delta, gamma):
                return -np.log10(10.**(-alp*x)+1)+delta*(np.log10(1.+np.exp(x)))**gamma/( 1. + np.exp(10.**-x) )
        
            def SHMRbeh(self, xmh,zz):
                nu, alp, delta, gamma, eps, xM1 = self.evolved_factors(zz)
                xmstel=eps+xM1+self.fbeh12(xmh-xM1, nu, alp, delta, gamma)-self.fbeh12(0.0, nu, alp, delta, gamma);
                return xmstel;
aa = Behroozi_2013a()
Mhalo = np.linspace(9.0, 15.0, 100)
#..................interpolation of sfe.dat file to find sfe=f(1+z,log M_hallo) ......
d_sfe=np.loadtxt('/home/iucaa/Downloads/sfh_z0_z8/release-sfh_z0_z8_052913/sfe/sfe.dat')
d_sfe0=d_sfe[:,0]
d_sfe1=d_sfe[:,1]
d_sfe2=d_sfe[:,2]
v=np.vstack((d_sfe0,d_sfe1)).T
sfe =interpolate.NearestNDInterpolator(v, d_sfe2)        
     


#----defining closest function to find redshift closer to redshift at which you required (zr) the plot--------

def closest(lst, K): 
      lst = np.asarray(lst) 
      idx = (np.abs(lst - K)).argmin() 
      return lst[idx]


#----reading tree files-------
drtry="/home/iucaa/Downloads/Stud15"
int_f=0  #initial file number


m_h=[]     #creates empty-list to save hallo mass at required z
m_s=[]     #creates empty-list to save stellar mass at required z

for FN in os.listdir(drtry):
    int_f=int_f+1
#-------to choose files only starts with name 'tree_'--------    
    if FN.startswith("tree_"):
        f=os.path.join(drtry,FN)
        d= np.loadtxt(f)        # d is the whole data set of one tree file

        print(int_f,FN)
        


#-------cutting data points from one tree file with z>5--and creating an array for redshift-----
        a=d[:,0]
        z_t=[]
        for i in range (len(d[:,0])):
            if ((1/a[i])-1<=5):
                z=(1/a[i])-1
        #        print(z)
                z_t.append(z)
                last_scale=i # --last row where z=5---
            else:
                break
        z_t=np.array(z_t,dtype=float)  #--------converting list to array-------
#-------storing data with z<=5-----   
        
        d1=d[0:last_scale+1]      #---d1 contains  data (tree hallos) with z<=5------
        
        scale_t=d1[:,0]           #----scale_t is scale factor.....
        num_prog=d1[:,1]          
        pid=d1[:,2]
        mvir=d1[:,3]              #----mvir is the given virial mass...
        
#----creating dictionary for a hallo containing all it's data...........       
        for i in range(len(scale_t)):
            f="h%d = {}"%i           # ----creats a dictionary (e.g 'h1' is a dictionary for first hallo)
            exec(f)               
            exec("h%d['scale'] = scale_t[i] "%i)   #--keys have there defined meaing.----
            exec("h%d['np'] = num_prog[i] "%i)
            exec("h%d['mvir'] = mvir[i] "%i)
            exec("h%d['z']=z_t[i]"%i)
        #print(h4722)
        zm=max(z_t)                #zm=maximum redshift
        
#----creating tree --- --
        #           -- - 
        #             -
#----assigning 0 no. of progenitor to  hallos of highest z----------   
        
        for i in range(len(z_t)):
            if(z_t[i]==zm):
                exec("h%d['np'] =0 "%i) #-- this will change the key value 'np'
               
#---------assigning progentors to hallos----------            
                
        k=1                            #--begins with 2nd hallo
        for i in range(len(scale_t)):
            n=[]                      #----a list to contain progenitor's hallo id (e.g ['h1','h2',..])
            exec("np1=h%d['np']"%i)
          
            for j in range(int(np1)):
                exec("p='h%d'"%k)       #--p is a string with named as corsponding hallo id (e.g p='h2')
                k=k+1
                n.append(p)
            exec("h%d['progenitors']=n"%i) #--we have a list for progenitors with elements as strings
                                           # (e.g h1['progenitors']=['h2','h3',....])
        
        

#.....Assigning stellar mass  to begining hallos and hallo with 0 no. of. progenitors.............
        for i in range(len(scale_t)-1,-1,-1):
            exec("np2=h%d['np']"%i)
            if(np2==0):
                mstelar = aa.SHMRbeh(np.log10(mvir[i]), z_t[i])
                exec("h%d['mstelar']=10**(mstelar)"%i)
#                exec("h%d['micl']=0"%i)
                
 #...................Finding ' dt'..................
        H_0=73.5e-12    # H_0  is Huble const. in unit of 1/years
 #--taking, (Omega_M=0.321),Omega_radiation(1e-4),Omega_darkenergy(0.678)------
        def dt(a,da):
            return da/(a*H_0*np.sqrt(0.321*a**(-3)+0.678+1e-4*a**(-4)))      
        
     #...........assigning dm/dt..................................
        for i in range(len(scale_t)-1,-1,-1):
            exec("np3=h%d['np']"%i)
            
    #------for 0 number of progenitors dm/dt=0---
            if(np3==0):
                exec("h%d['dm_dt'] =0 "%i)
    # ---finding dm_hallo for hallos with  non 0 np..........
            else:
                mh=np.zeros(int(np3))
                for j in range (int(np3)):
                    exec("g=h%d['progenitors'][j]"%i)
                    exec("m="+g)
                    mh[j]=m['mvir']
                mh_max=max(mh)
                scale_p=m['scale']
                exec("m_hallo=h%d['mvir']"%i)
                dm_h=m_hallo-mh_max
   #------finding dm/dt for dm_h >0 and dm/dt=0 for other dm_h 
                if (dm_h>0):
                    exec("scale_d=h%d['scale']"%i)
                    da=scale_d-scale_p  #  da is scale gap in between hallo scale and parent scale
                    dt_p=dt(scale_d,da) # dt is time gap in between hallo scale and parent scale
                    dmh_dt=dm_h/abs(dt_p) # determines dm_h/dt
                    sfe_d=sfe(1/scale_d,np.log10(mh_max)) # sfe of parent hallo with maximum virial mass in log form
                    exec("h%d['dm_dt']=dmh_dt*(10**sfe_d)*0.14"%i) #dm/dt=dm_h/dt. sfe . fb(=0.14)-----
                else:
                    exec("h%d['dm_dt']=0"%i)
        
        
   #-----assigning stellar mass to hallos with implimenting rules proposed by Becker...........        
        for i in range(len(scale_t)-1,-1,-1):
            
            n=int(locals()["h"+str(i)]['np'])  #--locals convert string 'hi' to hi{} dictionary-,n is np of ith hallo--
            a1=locals()["h"+str(i)]['scale']   #--a1 is scale of ith hallo 
    
            if n!=0:                           #for n=0 stellar mass previously assigned
                x=locals()["h"+str(i)]['progenitors'] #.x is the list containg  all progintiors of ith hallo in string format..

           #for n=1, h[M_stellar]=parent[M_stellar]+parent[dm/dt]*dt(time gap between hallo and parent) 
          
	        #---dm/dt now modified as dm/dt=dm/dt+log_normal (10^(-12)*[M_stellar],0.25)-------

                if n==1:               
                    j=locals()[x[0]]     # j is the parent dictionary
            
                    a2=j['scale']
                    da=(a1-a2)            #dt=f(a1,da) is time gap between parent -hallo and required hallo
                    m_star=j['mstelar']+(j['dm_dt']+(10**np.random.normal(np.log10(10**(-12)*j['mstelar']),0.25)))*dt(a1,da)
                    
                    locals()["h"+str(i)]['mstelar']=m_star
         # for np=2, stellar mass is added up along with newly formed stellar mass(i.e, dt*(parent['dm/dt']))
                if n==2:
                    j1=locals()[x[0]]
                    j2=locals()[x[1]]
                    a2=j1['scale']
                    da=(a1-a2)
                    m_star=j1['mstelar']+(j1['dm_dt']+10**np.random.normal(np.log10(10**(-12)*j1['mstelar']),0.25))*dt(a1,da)
                    m_star=m_star+j2['mstelar']+(j2['dm_dt']+10**np.random.normal(np.log10(10**(-12)*j2['mstelar']),0.25))*dt(a1,da)
                    locals()["h"+str(i)]['mstelar']=m_star

	# for np>2, stellar mass of massive progenitor along with newly formed stellar mass of massive progenitor only...
	
                if n>2:
               # finding massive progenitor id...
                    m_p=np.zeros(n)  
                    for k in range (n):
                        j=locals()[x[k]]
                        m_p[k]=j['mvir'] # m_p is an array, which contains mass of all progenitors..
                    max_mass=max(m_p) 
                    for k in range (n):
                        if (m_p[k]==max_mass):
                            n_max=k
                        else:
                            continue
                    j_max=locals()[x[n_max]] #j_max is the massive progenitor id..
                    a2=j_max['scale']
                    da=(a1-a2)
                    m_star=j_max['mstelar']+(j_max['dm_dt']+(10**np.random.normal(np.log10(10**(-12)*j_max['mstelar']),0.25)))*dt(a1,da)
                    
                    locals()["h"+str(i)]['mstelar']=m_star  #  you are done.stellar mass is assigned ..
#..to store the M_stellar and M_hallo for hallos with redshift closest to required redshift(zr).......

        zx=closest(z_t,zr)               # zx is the closest redshift to required z(zr) on current data file..

 
        for i in range(len(scale_t)):
            z_1=locals()["h"+str(i)]['z']
            if z_1==zx :
                m_h.append(locals()["h"+str(i)]['mvir'])
                m_s.append(locals()["h"+str(i)]['mstelar'])
                
m_h=np.array(m_h,dtype=float)    # converts list to array 
m_s=np.array(m_s,dtype=float)    
fig=plt.figure()
ax =fig.add_subplot(111)
Mhalo = np.linspace(9.0, 15.0, 100)       
Mstel = aa.SHMRbeh(Mhalo,zx)    # M_stellar from Behroozi function at closet redshift 'zx'
ax.plot(Mhalo, Mstel, label="z=%.3f" % zx)
ax.scatter(np.log10(m_h),np.log10(m_s),s=1)
ax.set_xlim(8.0, 16.0)
ax.set_ylim(7.0, 12.0)
ax.set_xlabel(r"$\log_{10} M_{\rm halo}$ ($h_{70}^{-1}M_\odot$)")
ax.set_ylabel(r"$\log_{10} M_*$ ($h_{70}^{-2}M_\odot$)")
ax.legend()
plt.show()
fig.savefig("z_%.1f.png"%zr)
end = it.time()
print(end - start)
