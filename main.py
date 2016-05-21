import matplotlib.pyplot as plt
from BlackSmith import BlackSmith


def main():
    ## Example execution
    nb_vertices = 1000
    mean_degs = [10,25,50] #mean degress for Gnp's
    nb_colors = 30
    CoolingPeriod = 20 #TODO: optimize, try different parameters.This parameter is really the unit for the x axis in the plot below
    T_guess = [30,40,50] #This parameter is very sensitive, if increase can augment running time a lot, but also tends to make the initial jump shorter

    #itertools.repeat(T_guess,len(mean_degs)
    ct_pairs = zip(mean_degs,T_guess)
    bsm = BlackSmith(q=nb_colors,n=nb_vertices,T_Period=CoolingPeriod,workline=ct_pairs,desc_thresh=0.8,schedule='Exp')

    for c,T,E,X,trace in bsm.work():
        plt.plot(trace,label='c=%d - T=%.2f'%(c,T))
    plt.legend(loc='best')
    plt.title('G(%d,c/%d) for a %d-coloring'%(nb_vertices,nb_vertices,nb_colors))

if __name__ == "__main__":
    main()
