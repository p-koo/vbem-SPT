function N = RandomTrackLength(numSim,meanN,minN,maxN)

rate = (maxN-minN)/(meanN-minN);
N = (maxN-minN)*exp(-rand(numSim,1)*rate)+minN;