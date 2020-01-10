def linearplot(x,y,varname):
    #do linear regression ,get the parameter and plot
    import matplotlib.pyplot as plt

    import numpy as np

    from sklearn import linear_model
    from sklearn.metrics import mean_squared_error, r2_score
    # x=np.linspace(1,10,10)
    # y=np.power(x,3)
    # varname='test'
    x=np.array(x).reshape(-1,1)
    y = np.array(y).reshape(-1, 1)
    regr=linear_model.LinearRegression()
    regr.fit(x,y)
    p=regr.predict(x)
    plt.scatter(x,y,c='b')
    minx=np.min(x)
    miny=np.min(y)
    maxx=np.max(x)
    maxy=np.max(y)
    minp=np.min(p)
    maxp=np.max(p)
    mmin=min(minx,miny,minp)
    mmax=max(maxx,maxy,maxp)
    r2=float('%.2f'%r2_score(y,p))
    rmse=float('%.2f'%np.sqrt(mean_squared_error(y,p)))
    n=len(x)

    plt.plot([mmin,mmax],[mmin,mmax],'k:')#1:1line
    plt.plot([minx,maxx],[minp,maxp],'r-')#regressionline

    plt.title(varname)
    plt.xlabel('QAAv6 derived/$m^{-1}$')
    plt.ylabel('in-situ/$m^{-1}$')
    #plt.xticks([])
    a=float('%.2f'%regr.coef_)
    b=float('%.2f'%regr.intercept_)
    plt.text(mmax-2,mmin,"$R^{2}=$"+str(r2)+"\n"+"RMSE="+str(rmse)+"\n"+"y="+str(a)+"x+"+str(b)+"\n"+"n="+str(n))




    plt.yscale('symlog')
    plt.xscale('symlog')

    plt.plot()
    plt.xticks([-0.5, 0.01, 0.5, 1,5])
    plt.yticks([-0.5, 0.01, 0.5, 1,5])
    plt.show()
    return [a,b,rmse,r2]
