#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:25:50 2018

@author: mlpoterek
"""

import pandas
import scipy
import scipy.integrate as spint
from plotnine import *

# Challenge - Lotka Volterra
# Assumption 1: only limiting factor for prey population growth is the predator
# Assumption 2: the only limiting factor for predator pop growth is the availibility of prey
def RMSim(y,t0,b,e,a,s,w,d):
    H=y[0]
    P=y[1]
    dHdt=H*b*(1-a*H) - w*(H/(d+H))*P
    dPdt=P*e*w*(H/(d+H)) - P*s
    #print(dHdt)
    
    return [dHdt,dPdt]
'''
# Case 1
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
simDF=pandas.DataFrame({"t":times,"prey":sim[:,0],"predator":sim[:,1]})
p=ggplot(simDF,aes(x="t",y="prey"))+geom_line(size=1,color='green')+geom_line(simDF,aes(x="t",y="predator"),color="red",size=1)
p=p+labs(title="Rosenzweig-MacArthur Model \n Suggested Parameters",x="Time",y="Population")+theme_classic()
p=p+scale_x_continuous(breaks=(0,100,200,300,400,500))+scale_y_continuous(breaks=(0,200,400,600,800,1000))
print(p)
'''


# Paradox of Enrichment
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
simDF=pandas.DataFrame({"t":times,"prey":sim[:,0],"predator":sim[:,1]})
q=ggplot(simDF,aes(x="t",y="prey"))+geom_line(size=1,color='green')+geom_line(simDF,aes(x="t",y="predator"),color="red",size=1)
q=q+labs(title="Rosenzweig-MacArthur Model \n alpha=0.001",x="Time",y="Population")+theme_classic()
q=q+scale_x_continuous(breaks=(0,100,200,300,400,500))+scale_y_continuous(breaks=(0,200,400,600,800,1000,1200,1400,1600,1800,2000))
print(q)
#q.save('RM_a=0.001.png')