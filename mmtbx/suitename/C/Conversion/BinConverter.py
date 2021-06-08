from pyparsing import nestedExpr
import suitenamedefs
from suitenamedefs import Bin, Cluster, SatelliteInfo

# paste the C initializer block for one bin in the triple quotes below
# then run.
# output is the python data block, as nested tuples
txt = '''
      {/*bin 13 incomplete angles*/
         "inc ",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "__", 0, "incompl", "white      ", "inc", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
'''

def isSequence(obj):
    return type(obj) in [list, tuple]


def commentize(list1):
    stringy = ""
    commented = False
    list2 = []
    for x in list1:
        if isSequence(x):
            list2.append(commentize (x))
        elif commented:
            if x[-2:] == "*/":
                stringy += x[:-2]
                commented = False
                list2.append(stringy)
                stringy = ""
            else:
                stringy += x + ' '
        elif x[:2] == "/*":
            commented = True
            stringy = x[2:] + ' '
        elif x == ",":
            pass
        else:
            list2.append(x)
    return list2


def convert(list1):
    list2 = list1[0]
    name = list2[0]
    i = int(name[4:6])
    print(f"bin{i}data = ({list2[1]},")
    list3 = list2[2]
    for i, list4 in enumerate(list3[1:], 1):
        ordinal = list4[0].strip()
        if ordinal == "N":
            break
        print("    (", end="")
        #list4[2] = (list4[2][1] != "0") # 0->False, 1->True
        outputSequence(list4[0:2])
        outputSequence(list4[3:6])
        print("\n        ", end="")
        outputNumbers(list4[7:])
        print("),")
    print(")\n\n")


# This alternative calls the cluster constructor on each line
def convertv2(list1, i):
    list2 = list1[0]
    name = list2[0]
    print(f"bins[{i}] = Bin({list2[1]})")
    print(f"bins[{i}].cluster = (")
    list3 = list2[2]
    for i, list4 in enumerate(list3[1:], 1):
        if list4[0][1] == "N":
            break
        print("    Cluster(", end="")
        list4[2] = (list4[2][1] != "0") # 0->False, 1->True
        outputSequence(list4[0:6])
        print("\n        ", end="")
        outputNumbers(list4[7:])
        print("),")
    print(")\n\n")


def outputSequence(a):
    for x in a:
        print(x, ", ", sep="", end="")


def outputNumbers(a):
    print("(", end="")
    if len(a) == 1:
        a = 9 * ('0',)
    for i, x in enumerate(a):
        if x[-1] == ",":
            x = x[:-1]
        elif x[0]==",":
            x = x[1:]
        if i==0:
            print(f"{float(x):07.3f}",end="")
        else:
            print(f", {float(x):07.3f}", end="")
    print(")", end="")

bin8data = ("23 t",
    ( 0 , "!!", "outlier", "white      ", "out", 
        (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
    ( 1 , "2h", "certain", "sea        ", "ord", 
        (180.0,  147.782,  260.712,  290.424,  296.2,  177.282,  175.594,  86.565,  180.0)),
    ( 2 , "4n", "certain", "peach      ", "ord",
        (180.0,  143.722,  227.256,  203.789,  73.856,  216.733,  194.444,  80.911,  180.0)),
    ( 3 , "0i", "certain", "lilactint  ", "sat",
        (180.0,  148.717,  274.683,  100.283,  80.6,  248.133,  181.817,  82.6,  180.0)),
    ( 4 , "6n", "certain", "lilac      ", "dom",
        (180.0,  150.311,  268.383,  84.972,  63.811,  191.483,  176.644,  85.6,  180.0)),
    ( 5 , "6j", "certain", "purple     ", "sat",
        (180.0,  141.633,  244.1,  66.056,  71.667,  122.167,  182.2,  83.622,  180.0)),
)

def buildBin(data):
    name = data[0]
    clusters = []
    for item in data[1:]:
        c = suitenamedefs.Cluster(*item)
        clusters.append(c)
    bin = Bin(name, clusters)
    return bin


out = nestedExpr('{','}').parseString(txt).asList()
out2 = commentize(out)
print(out2)
convert(out2)

# bin8 data is a sample of the output of convert
def formatIntegers(list):
    texts = [f'{n:2}' for n in list]
    text = ', '.join(texts)
    return text


def satelliteConvert(data):
    for k in range(0, 180, 20):
        name = data[k]
        satWidths = data[k+1:k+10]
        domWidths = data[k+10:k+19]
        doma = data[19]
        line = f'("{name}", ({formatIntegers(satWidths)}), ({formatIntegers(domWidths)}), "{doma}"),'
        print(line)


sats = (
    #   sat 0 1 s2 s3 4 s5 6 7 8 0 1 d2 d3 4 d5 6 7 8   dom
       "1m",0,0, 0, 0,0,32,0,0,0,0,0, 0, 0,0,64,0,0,0 ,"1a",
       "1L",0,0,18, 0,0,18,0,0,0,0,0,70, 0,0,70,0,0,0 ,"1a",
       "&a",0,0,20,20,0, 0,0,0,0,0,0,60,60,0, 0,0,0,0 ,"1a",
       "1f",0,0, 0, 0,0,47,0,0,0,0,0, 0, 0,0,65,0,0,0 ,"1c",
       "1[",0,0, 0, 0,0,34,0,0,0,0,0, 0, 0,0,56,0,0,0 ,"1b",
       "4a",0,0,40,40,0, 0,0,0,0,0,0,50,50,0, 0,0,0,0 ,"0a",
       "#a",0,0,26,26,0, 0,0,0,0,0,0,36,36,0, 0,0,0,0 ,"0a", 
       "0i",0,0, 0, 0,0,60,0,0,0,0,0, 0, 0,0,60,0,0,0 ,"6n",
       "6j",0,0, 0, 0,0,60,0,0,0,0,0, 0, 0,0,60,0,0,0 ,"6n"
)
# print(satelliteConvert(sats))

