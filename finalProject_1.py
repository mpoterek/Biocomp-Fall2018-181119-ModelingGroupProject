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

# Rosenzweig-MacArthur Model
def RMSim(y,t0,b,e,a,s,w,d):
    H=y[0]
    P=y[1]
    dHdt=H*b*(1-a*H) - w*(H/(d+H))*P
    dPdt=P*e*w*(H/(d+H)) - P*s
    
    return [dHdt,dPdt]

# Case 1
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
simDF=pandas.DataFrame({"t":times,"prey":sim[:,0],"predator":sim[:,1]})
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
p=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
p=p+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
p.save('RMstandard.png')


# Case 2: altered the prey birth rate "b"
times=range(1,500)
y0=[500.,120.]
parameters=(.3,0.07,0.001,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
q=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
q=q+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
q.save('RMdecreasedB.png')
#increasing the prey birth rate "b" did not change the final prey population,
#likely because of the system's carrying capacity, but it did increase the 
#predator population
#decreasing the prey birth rate "b" decreased the maximum population of the 
# herbivore and the max population of the predator by similar proportions


# Case 3: altered the predator attack rate "a"
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.0005,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
r=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
r=r+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
r.save('RMdecreasedA')
#decreasing the predator attack rate "a" caused oscillations in both populations;
#the prey population oscillated between a number in excess of the system's likely
#carrying capacity and ~0, while the predator population experienced much smaller
#oscillations
#increasing the predator attack rate "a" decreased the final prey and predator 
#populations


# Case 4: altered the conversion efficiency "e" of prey to predators
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.1,0.001,0.2,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
s=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
s=s+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
s.save('RMincreasedE')
#Increasing the conversion efficiency caused oscillations in both populations;
#the peak of the oscillations in both populations was above previous final values
#and the low point in both populations were close to zero
#Decreasing the conversion efficiency increased the final prey population and
#decreased the final predator population


# Case 5: altered the predator death rate "s"
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.1,5,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
t=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
t=t+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
t.save('RMdecreasedS')
#Decreasing the death rate caused oscillations in both populations; the peak
#of the oscillations in both populations was above previous final values and the
#low point in both populations were close to zero
#Increasing the predator death rate increased the final prey population and
#decreased the final predator population


# Case 6: altered the prey self-limitation factor "w"
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.2,7,400)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
u=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
u=u+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
u.save('RMincreasedW')
#Increasing w caused oscillations in both populations; the peak of the 
#oscillations in both populations was slightly above previous final values and 
#the low point in both populations were close to zero
#Decreasing w increased the final prey population and decreased the final 
#predator population


# Case 7: altered the predator saturation rate "d"
times=range(1,500)
y0=[500.,120.]
parameters=(0.8,0.07,0.001,0.2,5,700)
sim=spint.odeint(func=RMSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
v=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
v=v+labs(title="Rosenzweig-MacArthur Model",x="Time",y="Population")+theme_classic()
v.save('RMincreasedD')
#Increasing d increased the final prey population and decreased the final 
#predator population 
#Decreasing d caused oscillations in both populations; the peak of the 
#oscillations in the prey population was above the previous final value and the 
#low point in both populations were close to zero

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