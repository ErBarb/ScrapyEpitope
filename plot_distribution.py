import matplotlib.pyplot as plt
fp = open('/content/drive/My Drive/Colab Notebooks/Data/HPLG/epitope_table_export_1654247913.csv','r') # Link to file containing epitopes data
start_pos = []
end_pos = []
for ln in fp:
  # print(ln)
  sentence = ln.replace('\n','')
  sentence = sentence.split(',')

# Ignoring discontinuous epitopes and start and end position of some epitopes are blanks, so removing them
  if (sentence[1] == '"Linear peptide"') and (sentence[5] != '""') and (sentence[6] != '""'):
    sentence[5] = sentence[5].replace('\"','')
    sentence[6] = sentence[6].replace('\"','')

  # print(sentence)
    start_pos.append(sentence[5])
    end_pos.append(sentence[6])

# start_pos.pop()
# end_pos.pop()
# del start_pos[0]
# del start_pos[0]
# del end_pos[0]
# del end_pos[0]

print(len(start_pos))
print(len(end_pos))
print(start_pos)
print(end_pos)

count = [0]*1273
# for i in range(0,1272):
#   count[i].append(0)
print(len(count))
print(count)

for i in range(0,len(start_pos)):
  start_pos[i] = int(start_pos[i])
  end_pos[i] = int(end_pos[i])
  # print(i)
  # print(start_pos[i])
  # print(end_pos[i])
  for j in range(start_pos[i],end_pos[i]):
    count[j] = count[j]+1
print(count)

x_axis = [i for i in range(1, 1274)]
# for i in range(1,1274):
#   x_axis.append(i)
x_axis = list(range(1,1274))
print(len(x_axis))
print(x_axis)

plt.plot(x_axis,count)
# plt.title('title name')
plt.xlabel('Spike protein (P0DTC2) residues')
plt.ylabel('#Epitopes')
# plt.show()
plt.savefig('/content/drive/My Drive/Colab Notebooks/Data/plot1.png', dpi=300)
