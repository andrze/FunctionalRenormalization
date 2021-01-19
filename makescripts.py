from os import system
from sys import argv

mathematica_path = "/opt/Wolfram/Mathematica/12.0.0/Executables/math"

def run(mode, regulator, num_conf, additional=""):

    
    if mode not in ("tricritical","Z4", "ON", "lowd", "regulator"):
        raise ValueError("Invalid task specification, type: tricritical, Z4 or ON")
    if regulator not in ("smooth", "exponential"):
        raise ValueError("Invalid regulator specification, type: smooth or exponential")
    
    home = "/home/2/ac357729/Documents/FunctionalRenormalization"
    
    
    for i in range(1,num_conf+1):
        filename = "%s%02i.sh" % (mode,i)
        file = open("%s/tasks/%s" %(home,filename), 'w')   
        
        file.write("cd %s\n" % home)
        file.write("%s -script Run.wl -regulator %s -task %s -conf %i %s\n" % (mathematica_path, regulator, mode, i, additional))
        file.close()
    
        command = 'qsub -l mem=8000mb,walltime=240:00:00 %s/tasks/%s' % (home, filename)
        print command
        system(command)

if len(argv) == 4:
    run(argv[1],argv[2],int(argv[3]))
elif len(argv) > 4:
    run(argv[1],argv[2],int(argv[3]),argv[4])


