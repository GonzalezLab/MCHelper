# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import pandas as pd
import sys
from statistics import mean


fig = plt.figure(figsize=(20, 10))

# Creating dataset
fasta_file = sys.argv[1]
output_path = sys.argv[2]
output_file = sys.argv[3]


data_1 = []
data_2 = []
data_3 = []
data_4 = []
for te in SeqIO.parse(fasta_file, "fasta"):
    if "#" in str(te.id):
        order = str(te.id).split("#")[1].split(" ")[0]
        if "LTR" in order.upper():
            data_1.append(len(str(te.seq)))
        elif "LINE" in order.upper() or "SINE" in order.upper():
            data_2.append(len(str(te.seq)))
        elif "HELITRON" in order.upper():
            data_4.append(len(str(te.seq)))
        elif "TIR" in order.upper() or "DNA" in order.upper() or "MITE" in order.upper():
            data_3.append(len(str(te.seq)))


mean1 = 0
mean2 = 0
mean3 = 0
mean4 = 0
if len(data_1) > 0:
    mean1 = mean(data_1)
if len(data_2) > 0:
    mean2 = mean(data_2)
if len(data_3) > 0:
    mean3 = mean(data_3)
if len(data_4) > 0:
    mean4 = mean(data_4)

print(str(mean1) + ";" + str(mean2) + ";" + str(mean3) + ";" + str(mean4))

dataframe = pd.DataFrame([data_1, data_2, data_3, data_4])
dataframe.to_csv(fasta_file+"_length.csv", sep=";")

# Creating plot
data = [data_1, data_2, data_3, data_4]
fig, ax = plt.subplots()
# Creating axes instance
bp = ax.boxplot(data)
ax.set_ylim(0, 25000)
# x-axis labels
ax.set_xticklabels(['LTR', 'LINE', 'TIR', 'Helitron'])
#plt.show()
plt.savefig(output_path+'/'+output_file)

