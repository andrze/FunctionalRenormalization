from os import system
from sys import argv

mathematica_path = "/opt/Wolfram/Mathematica/12.0.0/Executables/math"
time = {"ising" : 48,
        "Z4" : 120,
        "ON" : 240}

def run(mode, regulator, num_conf, additional="", start_index=1):

    tasks = ("tricritical","Z4", "ON", "lowd", "regulator", "ising", "ising2")
    if mode not in tasks:
        raise ValueError("Invalid task specification, type one of: " + ", ".join(tasks))
    if regulator not in ("smooth", "exponential", "litim"):
        raise ValueError("Invalid regulator specification, type: smooth or exponential")
    
    home = "/home/2/ac357729/Documents/FunctionalRenormalization"
    
    memory = 8000
    if mode=="ON":
        memory = 24000
    
    for i in range(start_index,num_conf+1):
        if num_conf < 1000:
            filename = "%s%03i.sh" % (mode,i)
        else:
            filename = "%s%04i.sh" % (mode,i)
        
        file = open("%s/tasks/%s" %(home,filename), 'w')   
        
        file.write("cd %s\n" % home)
        file.write("%s -script Run.wl -regulator %s -task %s -conf %i %s\n" % (mathematica_path, regulator, mode, i, additional))
        file.close()
    
    
        hours = time.get(mode, 120)
        command = 'qsub -l mem=%imb,walltime=%i:00:00 %s/tasks/%s' % (memory, hours, home, filename)
        print command
        system(command)

if len(argv) == 4:
    run(argv[1],argv[2],int(argv[3]))
elif len(argv) == 5:
    run(argv[1],argv[2],int(argv[3]),argv[4])
elif len(argv) > 5:
    run(argv[1],argv[2],int(argv[3]),argv[4], int(argv[5]))


