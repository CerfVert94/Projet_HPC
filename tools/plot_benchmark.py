#!/usr/bin/python3

import matplotlib, matplotlib.pyplot as plt, sys, os, numpy as np

# matplotlib.use('')

if __name__ == "__main__":
    # print(sys.argv[1])
    if not os.path.isfile(sys.argv[1]) :
        raise ValueError(sys.argv[1] + " doesn't exist.")
    file = open(sys.argv[1], "r")

    buffer = file.read()
    # stripped = buffer.replace("  ", " ")
    # while stripped != buffer :
    #     buffer = stripped
    #     stripped = buffer.replace("  ", " ")

    buffer = buffer.split('\n')
    labels = []

    buffer[0] = buffer[0].strip()
    # print(buffer[0])
    
    while len(buffer[0]) > 0 :
        labels.append(buffer[0].split(' ')[0])
        buffer[0] = buffer[0][len(labels[-1]):].strip()
        # print(labels[-1], ";")

    size = np.zeros(len(buffer) - 1, dtype=int)
    data = np.zeros([len(buffer) - 1, len(labels) - 1])

    for i in range(1, len(buffer) - 1):
        buffer[i] = buffer[i].strip()
        size_str = buffer[i].split(' ')[0]
        buffer[i] = buffer[i][len(size_str) :].strip()
        size[i] = int(size_str)
        # print(size[i], end = ' ')
        for j in range(0, len(labels) - 1):
            cycles_str = buffer[i].split(' ')[0]
            buffer[i] = buffer[i][len(cycles_str) :].strip()
            data[i, j] = float(cycles_str)
            # print(data[i, j], end = ' ')
        # print(end ='\n')
        # input()
    
    best_labels = []
    fig = plt.figure()
    for i in range(0, len(labels) - 1):
        if i >= 0 and np.average(data[ 1:, i]) < 2.5:
            plt.plot(size, data[ :, i])
            plt.xlabel("Size (N x N)")
            plt.ylabel("Cycles per Point")
            print(labels[1 + i],":")
            best_labels.append(labels[1 + i])
            print("\tAvg : ", np.average(data[ 1:, i]))
            print("\tMin : ", np.min(data[ 1:, i]))
            print("\tMax : ", np.max(data[ 1:, i]))
    plt.legend(best_labels)
    plt.show()
            
        # print(buffer[i])


    

        
    
    file.close()




