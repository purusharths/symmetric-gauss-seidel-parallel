import subprocess
import itertools
# stencil_size = [2,4,8,12,26,20,24,30]
# block_size = [2,5,10,20,25]
# stencil_size = [3,4,7,8,15,16,20,24,30,32]
# block_size =   [2,5,10,20,25]
k_ = [5,10,15]
n_ = [100,200,300]
iterations = 100
for k,n in zip(k_,n_):
    print("Started Stencil of size n={} and k= {} ... ".format(k, n))
    args = ("./sgs", str(n), str(k), str(iterations))
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    print(output)
    break
    # with open("benchmark_results/{}_stencil_{}_block".format(stencil, block),"w+") as f:
    #     f.write(output.decode('ascii'))
print("Done")

