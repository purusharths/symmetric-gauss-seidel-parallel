import pandas as pd
import numpy as np
import sys

from matplotlib import pyplot as plt
from matplotlib import style

filename = sys.argv[1]
# title = sys.argv[2]

style.use('seaborn-colorblind')
df = pd.read_csv(filename, comment='#')

def make_dfs(dframe, col_name):
    unique_exps = dframe[col_name].value_counts().to_frame().index.to_list()
    unique_dfs = [dframe[dframe[col_name] == exp] for exp in unique_exps]
    return unique_exps, unique_dfs


y = ['time', 'Gflops/s', 'GByte/s']
x = 'n'

fig, axes = plt.subplots(nrows=1, ncols=len(y), figsize=(20, 6))
unique_exps, unique_dfs = make_dfs(dframe=df, col_name='experiment')
x_ = [each_df[x].to_list() for each_df in unique_dfs][0]


for i in range(len(y)):
    y_ = [each_df[y[i]].to_list() for each_df in unique_dfs]
    axes[i].set_title(y[i])

    for (exp, each_y) in zip(unique_exps, y_):
        axes[i].plot(x_, each_y, '^-', label=exp)
        axes[i].set_xlabel(x)
    axes[i].legend()
    axes[i].grid(True)
    # axes[i].title()

plt.show()
