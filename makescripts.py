from os import system
from sys import argv

mathematica_path = "/opt/Wolfram/Mathematica/12.0.0/Executables/math"

def run(mode, regulator, num_conf):

    
    if mode not in ("4", "ON"):
        raise ValueError("Invalid model specification, type: 4 or ON")
    if regulator not in ("smooth", "exponential"):
        raise ValueError("Invalid regulator specification, type: smooth or exponential")
    
    home = "/home/2/ac357729/Documents/FunctionalRenormalization"
    
    label = "Z4" if mode == "4" else "ON"
    
    for i in range(1,num_conf+1):
        filename = "%s%02i.sh" % (label,i)
        file = open("%s/tasks/%s" %(home,filename), 'w')   
        
        file.write("cd %s\n" % home)
        file.write("%s -script Run.wl regulator %s Zp %s conf %i\n" % (mathematica_path, regulator, mode, i))
        file.close()
    
        command = 'qsub -l mem=8000mb,walltime=48:00:00 %s/tasks/%s' % (home, filename)
        print command
        system(command)

run(argv[1],argv[2],int(argv[3]))

