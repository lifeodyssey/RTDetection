def errorplot(x,y,varname):
    #for y,column is the value in one band,row is the value in one bands
    #x is band
    import numpy as np
    import matplotlib.pyplot as plt
    L=len(x)
    mean=np.ones(L)
    error=np.ones(L)
    for i in range(L):
        mean[i]=np.mean(y[:,i])
        error[i]=np.std(y[:,i])

    plt.errorbar(x, mean, yerr=error)
    tick=np.arange(400,700,25)
    plt.xticks(tick)
    plt.xlabel("wavelength")
    plt.ylabel(varname)
    plt.title(varname)
    return [mean,error]