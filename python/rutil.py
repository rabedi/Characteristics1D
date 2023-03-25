import os

def rrep(str, repPart):
    if (repPart in str):
        strOut = str.replace(repPart, "")
        return strOut, True
    strOut = str
    return strOut, False

def changeName(str, repParts):
    for repPart in repParts:
        strNew, b = rrep(str, repPart)
        if b:
            return strNew

    strNew = str
    return strNew
        

def rmkdir(folderName):
    if os.path.exists(folderName):
        return
    try:
        os.mkdir(folderName)
        #print(f"Directory '{folderName}' created successfully!")
    except OSError as error:
        print(f"Directory '{folderName}' creation failed. Error: {error}")
