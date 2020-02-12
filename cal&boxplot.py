import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
file='I:/Nagoya University/Project/HAB/match up/all.xlsx'
df=pd.read_excel(file,sheet_name='Sheet5')
df.loc[:,'N412']=df["Rrs_412"]/df["Rrs_555"]
df.loc[:,'N443']=df["Rrs_443"]/df["Rrs_555"]
df.loc[:,'N469']=df["Rrs_469"]/df["Rrs_555"]
df.loc[:,'N488']=df["Rrs_488"]/df["Rrs_555"]
df.loc[:,'N531']=df["Rrs_531"]/df["Rrs_555"]
df.loc[:,'N547']=df["Rrs_547"]/df["Rrs_555"]
df.loc[:,'N555']=df["Rrs_555"]/df["Rrs_555"]
df.loc[:,'N645']=df["Rrs_645"]/df["Rrs_555"]
df.loc[:,'N667']=df["Rrs_667"]/df["Rrs_555"]
df.loc[:,'N678']=df["Rrs_678"]/df["Rrs_555"]
df.loc[:,'slope412_443']=(df["Rrs_412"]-df["Rrs_443"])/(412-443)
df.loc[:,'slope443_488']=(df["Rrs_443"]-df["Rrs_488"])/(443-488)
df.loc[:,'slope488_555']=(df["Rrs_488"]-df["Rrs_555"])/(488-555)
df.loc[:,'slope555_645']=(df["Rrs_555"]-df["Rrs_645"])/(555-645)
df.loc[:,'slope645_667']=(df["Rrs_645"]-df["Rrs_667"])/(645-667)
df.loc[:,'slope667_678']=(df["Rrs_667"]-df["Rrs_678"])/(667-678)
df.loc[:,'SS443']=(df["Rrs_443"]-df["Rrs_412"]-(df["Rrs_469"]-df["Rrs_412"])*((443-412)/(469-412)))
df.loc[:,'SS469']=(df["Rrs_469"]-df["Rrs_443"]-(df["Rrs_488"]-df["Rrs_443"])*((469-443)/(488-443)))
df.loc[:,'SS488']=(df["Rrs_488"]-df["Rrs_469"]-(df["Rrs_531"]-df["Rrs_469"])*((488-469)/(531-469)))
df.loc[:,'SS531']=(df["Rrs_531"]-df["Rrs_488"]-(df["Rrs_547"]-df["Rrs_488"])*((531-488)/(547-488)))
df.loc[:,'SS645']=(df["Rrs_645"]-df["Rrs_555"]-(df["Rrs_667"]-df["Rrs_555"])*((645-555)/(667-555)))
df.loc[:,'SS667']=(df["Rrs_667"]-df["Rrs_645"]-(df["Rrs_678"]-df["Rrs_645"])*((667-645)/(678-645)))

wavelength=[412,443,469,488,531,547,555,645,667,678]
#df.set_index('water type',inplace=True)

Nchattonella=df.loc[df['water type']=='chattonella','N412':'N678']
Ndiatom=df.loc[df['water type']=='diatom','N412':'N678']
Ndiatom.reset_index(drop=True, inplace=True)
Nothers=df.loc[df['water type']=='others','N412':'N678']
Nothers.reset_index(drop=True, inplace=True)
Nturbid=df.loc[df['water type']=='turbid','N412':'N678']
Nturbid.reset_index(drop=True, inplace=True)
plt.figure(1)
for i in range(Nchattonella.shape[0]):
    plt.plot(wavelength,Nchattonella.loc[i],'r',label='chattonella')

for i in range(5):
    plt.plot(wavelength,Ndiatom.loc[i],'m',label='diatom')

for i in range(Nothers.shape[0]):
    plt.plot(wavelength,Nothers.loc[i],'g',label='others')

for i in range(Nothers.shape[0]):
    plt.plot(wavelength, Nturbid.loc[i], 'peru',label='turbid')
plt.hlines(1,412,678,'r')
plt.vlines(wavelength,0,1.2,'k')
plt.xticks(wavelength)
plt.show()


Rchattonella=df.loc[df['water type']=='chattonella','Rrs_412':'Rrs_678']
Rdiatom=df.loc[df['water type']=='diatom','Rrs_412':'Rrs_678']
Rdiatom.reset_index(drop=True, inplace=True)
Rothers=df.loc[df['water type']=='others','Rrs_412':'Rrs_678']
Rothers.reset_index(drop=True, inplace=True)
Rturbid=df.loc[df['water type']=='turbid','Rrs_412':'Rrs_678']
Rturbid.reset_index(drop=True, inplace=True)
plt.figure(2)
for i in range(Rchattonella.shape[0]):
    plt.plot(wavelength,Rchattonella.loc[i],'g',label='chattonella')

for i in range(5):
    plt.plot(wavelength,Rdiatom.loc[i],'g',label='diatom')

for i in range(Rothers.shape[0]):
    plt.plot(wavelength,Rothers.loc[i],'g',label='others')

for i in range(Rothers.shape[0]):
    plt.plot(wavelength, Rturbid.loc[i], 'peru',label='turbid')
#plt.hlines(0.7,412,678,'r')
plt.vlines(wavelength,0,0.02,'k')
plt.xticks(wavelength)
plt.show()


plt.figure(3)
Slopechattonella=df.loc[df['water type']=='chattonella','slope412_443':'slope667_678']
Slopediatom=df.loc[df['water type']=='diatom','slope412_443':'slope667_678']
#Rdiatom.reset_index(drop=True, inplace=True)
# Rothers=df.loc[df['water type']=='others','Slope412_443':'Slope667_678']
# Rothers.reset_index(drop=True, inplace=True)
# Rturbid=df.loc[df['water type']=='turbid','Slope412_443':'Slope667_678']
# Rturbid.reset_index(drop=True, inplace=True)
Slopechattonella.boxplot()

Slopediatom.plot.box(color=dict(boxes='m',whiskers='m', medians='m', caps='m'))

plt.figure(4)
SSchattonella=df.loc[df['water type']=='chattonella','SS443':'SS667']
SSdiatom=df.loc[df['water type']=='diatom','SS443':'SS667']
SSchattonella.plot.box(color=dict(boxes='r',whiskers='r', medians='r', caps='r'))

SSdiatom.plot.box(color=dict(boxes='m',whiskers='m', medians='m', caps='m'))

plt.show()