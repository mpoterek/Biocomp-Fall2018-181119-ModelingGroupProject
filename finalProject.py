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
def LVSim(y,t0,b,a,e,s):
    H=y[0]
    P=y[1]
    dHdt=H*b - H*a*P
    dPdt=P*e*a*H - P*s
    #print(dHdt)
    
    return [dHdt,dPdt]

# Case 1
times=range(1,100)
y0=[25.,5.]
parameters=(0.5,0.02,0.1,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
simDF=pandas.DataFrame({"t":times,"prey":sim[:,0],"predator":sim[:,1]})
p=ggplot(simDF,aes("t","prey"))+geom_line(size=1,color='green')+geom_line(simDF,aes("t","predator"),color="red",size=1)
p=p+labs(title="Lotka-Volterra Model",x="Time",y="Population")+theme_classic()
print(p)
