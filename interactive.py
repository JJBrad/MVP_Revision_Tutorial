import sys

def readArgs(args, params):
    """
    Function which translates command line arguments and acts on them.
    :param args: List of arguments to the command line.
    """
    i=0
    updates = {}
    while i < len(args):
        if args[i] in ["help", "Help", "-h", "H", "h", "?", "??"]:
            with open("README.md", "r") as readme:
                print(readme.read())
                exit()
        elif args[i] in ["-RS", "-rs"]: 
            try:
                updates["Seed"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised seed.")
                exit()
        elif args[i] in ["-N", "-n"]: 
            try:
                updates["noise"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised seed.")
                exit()
        elif args[i] in ["-s", "-S"]: 
            try:
                updates["sigma"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -S.")
                exit()
        elif args[i] in ["-dt"]: 
            try:
                updates["dt"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -dt.")
                exit()
        elif args[i] in ["-dx"]: 
            try:
                updates["dx"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -dx.")
                exit()
        elif args[i] in ["-x", "-X"]:
            try:
                updates["X Dimension"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -x.")
                exit()
        elif args[i] in ["-T", "-t"]:
            try:
                updates["tMax"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -N.")
                exit()
        elif args[i] in ["-y", "-Y"]:
            try:
                updates["Y Dimension"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -y.")
                exit()
        elif args[i] in ["-D", "-d"]:
            try:
                updates["D"] = float(args[i+1])
                i += 2
            except:
                print("Error with -D tag.")
                exit()
        elif args[i] in ["-i", "-I"]:
            try:
                updates["phi0"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -I")
                exit()  
        elif args[i] in ["-K", "-k"]:
            try:
                updates["kappa"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -k")
                exit()
        elif args[i] in ["-u", "-U"]:
            try:
                updates["Rate"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -U")
                exit()
        elif args[i] in ["-P", "-p"]:
            try:
                updates["Interval"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -P")
                exit()
        elif args[i] in ["-v", "-V"]:
            try:
                updates["v0"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -v")
                exit()
        else:
            print("Key {} not recognised. Ignoring.".format(args[i]))
            i += 2
    params.update(updates)
        
        
                
            
            
