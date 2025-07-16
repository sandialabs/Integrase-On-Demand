
import os,sys,shutil,subprocess,time

start=time.time()
inputs = sys.argv               # get arguments and initialize defaults
defaults={"-yl":"10",
          "-sl":"16",
          "-threads" : "1",
          "-tax" : "False",
          "-search" : "False",
          "-seq" : "None",
          "-n": "500",
          "-logfile": "iod.log"
          }

############################################################################################################
######################################### usage & options ##################################################

usage="""python main.py -in <path to genome.fasta or GCA> -out <path to desired output directory> -tax | -search\nTo see optional arguments, use -h"""

options= """OPTIONS\n-out <path to desired output directory>, default='./Integrase_on_Demand_Output/'
\n-yl <integer> length of matching flanks for tyrosine integrase attB sites, default=10
\n-sl length of matching flanks for serine integrase attB sites, default=16
\n-threads number of threads to use for blastn, default=1
\n-seq attB sequence or path to fasta containing sequences
\n-n max number of islands to allow in taxonomic mode, default=5000."""

modes="-tax (default)| -search"

mode_description="""MODES {modes}\n1. Taxonomic mode assigns a species to the input genome and collects
attB sequences from TIGER/Islander predicted islands from close relatives. Use case: First attempt search for bacterium; faster than
searching through all islands. If this mode returns unsatisfactory results, consider -search mode.\n2.Search mode takes a list of 
specified attB sequences. Users can use this mode in three ways: 1. provide sequences as comma seperated list in standard input, 
2. provide a path to a file containing each sequence on a newline or in FASTA format, or 3. no sequences specified so the whole database
will be used"""

###########################################################################################################
########################################input parsing#######################################################

if any((i in ('-h', '-help','--help') for i in inputs)):     #handles user --help or -h
    print(f"\nUsage:\n{usage}\n{inputs}\n{options}\n")
    sys.exit()
if not  '-in' in inputs:        #handles if user fails to provide FASTA/GCA
    print(f"Usage error:\n" + usage)
    print("Must specify path to genome fasta file with flag -in <FilePath>")
    sys.exit()
if not '-tax' in inputs and not '-search' in inputs:
    inputs+=["-tax"]
    #print("***\nMUST SPECIFY A MODE\n1. TAXONOMIC -tax\n- OR -\n2. SEARCH -search [-seq <comma separated list of attBs, or path to FASTA>]\n***")
    #print(mode_description)
    #sys.exit()

genomelength=0
fasta_path = inputs[inputs.index('-in') + 1]
if not os.path.exists(fasta_path):
    print(f"Could not find FASTA file in path: {fasta_path}\n{usage}")
    sys.exit()
with open(fasta_path,'r') as f:
    for i in f:
        if not i.startswith('>'):
            genomelength+=len(i)                      #check that only NGS bases in fasta
                    #check if output path was given and ensure proper format for later calls
if "-out" in inputs:
    i = inputs.index('-out')
    output_path= inputs[i+1]
    if output_path[-1]!= '/':
        output_path+='/'
    if not os.path.exists(output_path):
        try:
            os.mkdir(output_path) #make output directory
        except Exception as e:
            print(f"Exiting. Unable to create output directory due to the following system error: {e}")
            sys.exit()
elif not '-out' in inputs:
    output_path = './Integrase_On_Demand_Output/'
    if not os.path.exists(output_path):
        try:
            os.mkdir(output_path)
        except Exception as e:
            print(f'Exiting: Unable to write to current directory due to the following system error: {e}\nTry specifying an output directory with -out <DirectoryName>')
            sys.exit() 

def HandleOptionalArgs(flag): #process optional values or set variable to defaults
    if not flag.startswith('-'):
        flag=f"-{flag}"
    if flag in inputs:
        if flag=="-tax":
            defaults[flag]="True"
            return
        elif flag=='-search':
            i=inputs.index(flag)
            if flag==inputs[-1] or inputs[i+1] in defaults: # this may give an error if -search is last arg
                defaults[flag]="True"
                return
            defaults["-seq"]=inputs[i+1]
            defaults[flag]="True"
        i=inputs.index(flag)
        val=inputs[i+1]  #captures values
        defaults[flag]=str(val)
    return

for i in defaults:
    HandleOptionalArgs(i)
runiod_path=sys.argv[0]
script_path=runiod_path.split("bin")[0]
if not os.path.exists(f"{script_path}bin/iod.py"):
   print(f"Exiting. iod.py was not found in {script_path}. Please move iod.py to the same directory as runiod.py({script_path}bin).")
   sys.exit()
if output_path.strip('./') in defaults["-logfile"].split('/'):
    full_log_path=defaults["-logfile"]
else:
    full_log_path=os.path.join(output_path,defaults["-logfile"])
del defaults["-logfile"]
if os.path.exists(full_log_path):
    os.remove(full_log_path)
os.system(f"echo 'Running Command: python {' '.join(sys.argv)}' > {full_log_path}")
os.system(f"echo 'Running Command: python -u {script_path}bin/iod.py {fasta_path} {output_path} {script_path}lib/ {' '.join(list(defaults.values()))} &>> {full_log_path}' >> {full_log_path}")
os.system(f"python -u {script_path}bin/iod.py {fasta_path} {output_path} {script_path}lib/ {' '.join(list(defaults.values()))} &>> {full_log_path}")
runtime=round(time.time()-start, 2)
os.system(f"echo 'Runtime={runtime}, Genome size={genomelength}, Rate={genomelength/runtime} bp/s' >> {full_log_path}")

