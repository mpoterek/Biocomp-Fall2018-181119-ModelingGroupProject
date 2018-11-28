# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 11:15:21 2018

Look at Exercise 10 part 2 to find a basis for the Lotka-Volterra model

@author: Patricia
"""

# load packages
import pandas
import scipy
import scipy.integrate as spint
from plotnine import *

# Challenge - Lotka Volterra
# Assumption 1: only limiting factor for prey population growth is the predator
# Assumption 2: the only limiting factor for predator pop growth is the availibility of prey
def LVSim(y,t0,H,P,b,a,e,s):
    N1=y[0]
    N2=y[1]
    dHdt=H*(b - a*P)*N1
    dPdt=P*(e*a*H - s)*N2
    
    
    return [dHdt,dPdt]

# Case 1
times=range(1,50)
y0=[0.1,0.1]
parameters=(25,5,0.5,0.02,0.1,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
#simDF=pandas.DataFrame({"t":times,"prey":sim[:,0],"predator":sim[:,1]})
#print(ggplot(simDF,aes(x="t",y="prey"))+geom_line()+geom_line(simDF,aes(x="t",y="predator"),color="red")+theme_classic())

