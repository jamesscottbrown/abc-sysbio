


def evaluateDistance(distance,epsilon,t):

    accepted = False
    
    #only 1 epsilon sequences
    if(len(epsilon)==1):
        if(distance<epsilon[0][t] and distance>=0 ):
            accepted = True

    #for several epsilon sequences
    else:
        for i in range(len(epsilon)):
            if(distance[i]<epsilon[i][t] and distance[i]>=0 ):
                accepted = True
            else: 
                accepted = False
                break

    return accepted
