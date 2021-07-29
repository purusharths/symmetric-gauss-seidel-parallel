import subprocess
import itertools
# stencil_size = [2,4,8,12,26,20,24,30]
# block_size = [2,5,10,20,25]
stencil_size = [3,4,7,8,15,16,20,24,30,32]
block_size =   [2,5,10,20,25]
for stencil,block in list(itertools.product(stencil_size, block_size)):
    print("Started Stencil of size n={} and block= {} ... ".format(stencil, block), end="")
    args = ("./benchmark", "10000", str(stencil), str(block))
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    with open("benchmark_results/{}_stencil_{}_block.csv".format(stencil, block),"w+") as f:
        f.write(output.decode('ascii'))
    print("Done")

