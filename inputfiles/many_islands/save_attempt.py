import os
import shutil


attemptsFolder = os.listdir('attempts/')
n = 0
for s in ['table', 'regions', 'B_ext', 'log']:
    n = max(n, sum(s in i for i in attemptsFolder)) # This is the %06.d number which we will use now
print(n)

os.system('mumax3-convert -png -arrows 8 "many_islands_interaction.out/*.ovf"')
outputFolder = os.listdir('many_islands_interaction.out/')
for s in ['table.txt', 'regions000000.png', 'B_ext000000.png', 'log.txt']: 
    if s not in outputFolder:
        raise FileExistsError('Could not save the attempt, because "many_islands_interaction.out/%s" does not exist.' % s)


# Copy files to attempts folder, with an incremented index n
shutil.copyfile('many_islands_interaction.out/table.txt', 'attempts/table%06.d.txt' % n)
shutil.copyfile('many_islands_interaction.out/regions000000.png', 'attempts/regions%06.d.png' % n)
shutil.copyfile('many_islands_interaction.out/B_ext000000.png', 'attempts/B_ext%06.d.png' % n)
shutil.copyfile('many_islands_interaction.out/log.txt', 'attempts/log%06.d.txt' % n)

print('Successfully saved attempt as %06.d.' % n)