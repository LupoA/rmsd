import os
import h5py


T = 32


directory = r'mesons'
for filename in os.listdir(directory):
    if filename.endswith(".h5"):
        with open('PP.txt', 'a') as output:
            f = h5py.File(os.path.join(directory, filename), "r")
            #for key in f.keys():
            #    print(key) #Names of the groups in HDF5 file.
            #group = f[key]
            #for key in group.keys():
            #  print(key)
            for i in range(0,T):
                if i in range(0,T-1):
                    print( f['meson']['meson_0']['corr'][i][0], end="//", file=output)
                if i==T-1:
                    print( f['meson']['meson_0']['corr'][i][0], file=output)
        #rint(os.path.join(directory, filename))
    else:
        continue
print("Writing is over")