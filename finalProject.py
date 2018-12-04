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
import matplotlib.pyplot as plt
from plotnine import *
import numpy as np

# Challenge - Lotka Volterra
# Assumption 1: only limiting factor for prey population growth is the predator
# Assumption 2: the only limiting factor for predator pop growth is the availibility of prey
def LVSim(y,t0,b,a,e,s):
    H=y[0]
    P=y[1]
    dHdt=H*b - H*a*P
    dPdt=P*e*a*H - P*s
    
    return [dHdt,dPdt]

# Case 1: suggested parameters
times=range(1,100)
y0=[25.,5.]
parameters=(0.5,0.02,0.1,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
p=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
p=p+labs(title="Lotka-Volterra Model \n Suggested Parameters",x="Time",y="Population")+theme_classic()
p=p+scale_x_continuous(breaks=(0,10,20,30,40,50,60,70,80,90,100))+scale_y_continuous(breaks=(0,100,200,300,400,500))
print(p)
#p.save('LVstandard.png')




# Case 2: altered the prey birth rate "b"
times=range(1,100)
y0=[25.,5.]
parameters=(0.25,0.02,0.1,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
q=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
q=q+labs(title="Lotka-Volterra Model \n Decreased Herbivore Birth Rate",x="Time",y="Population")+theme_classic()
q=q+scale_x_continuous(breaks=(0,10,20,30,40,50,60,70,80,90,100))
print(q)
#q.save('LVdecreasedB.png')
# increasing the prey birth rate "b" increased the population of the herbivore
# as well as the population of the predator by similar proportions
# decreasing the prey birth rate "b" decreased the maximum population of the 
# herbivore and the max population of the predator by similar proportions


# Case 3: altered the predator attack rate "a"
times=range(1,100)
y0=[25.,5.]
parameters=(0.5,0.01,0.1,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
r=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
r=r+labs(title="Lotka-Volterra Model \n Decreased Predator Attack Rate",x="Time",y="Population")+theme_classic()
r=r+scale_x_continuous(breaks=(0,10,20,30,40,50,60,70,80,90,100))
print(r)
#r.save('LVdecreasedA')
# decreasing the predator attack rate "a" increased the maximum prey and predator populations
# but more importantly in the same time frame as the baseline, decreasing the 
# predator attack rate decreased the number of full cycles. 
# increasing the predator attack rate "a" decreased the maximum prey and predator populations
# and increased the number of cycles in the same time period


# Case 4: altered the conversion efficiency "e" of prey to predators
times=range(1,100)
y0=[25.,5.]
parameters=(0.5,0.02,0.05,0.2)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
s=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
s=s+labs(title="Lotka-Volterra Model \n Decreased Conversion Efficiency",x="Time",y="Population")+theme_classic()
s=s+scale_x_continuous(breaks=(0,10,20,30,40,50,60,70,80,90,100))
print(s)
#s.save('LVdecreasedE')
# Increasing the conversion efficiency did not seem to affect the max predator population
# but it lowered the prey population
# Decreasing the conversion efficiency also did not affect the max predator population
# but it raised the maximum prey population


# Case 5: altered the predator death rate "s"
times=range(1,100)
y0=[25.,5.]
parameters=(0.5,0.02,0.1,0.1)
sim=spint.odeint(func=LVSim,y0=y0,t=times,args=parameters)
modelOutputH=pandas.DataFrame({"t":times,"pop":list("H")*len(times), "density":sim[:,0]})
modelOutputP=pandas.DataFrame({"t":times,"pop":list("P")*len(times), "density":sim[:,1]})
modelOutput=pandas.concat([modelOutputH, modelOutputP])
t=ggplot(modelOutput,aes(x="t",y="density"))+geom_line(aes(color="pop"), size=1)+theme_classic()
t=t+labs(title="Lotka-Volterra Model \n Decreased Predator Death Rate",x="Time",y="Population")+theme_classic()
t=t+scale_x_continuous(breaks=(0,10,20,30,40,50,60,70,80,90,100))
print(t)
#t.save('LVdecreasedS')
# increasing the predator death rate did not affect the max predator population
# but it increased the maximum prey population as well as increased the number
# of cycles in the same time period
# decreasing the predator death rate again did not affect the max predator population
# but lwoered the max prey population as well as decreased the number of 
# cycles in the same time period