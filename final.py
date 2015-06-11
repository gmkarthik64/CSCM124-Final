import copy
import random
import itertools
import time
    
h = [None]*4
h[0] = []
h[1] = []
h[2] = []
h[3] = []
snplen = 100

reads = []

for o in range(snplen):
    i = random.randint(1,14)
    h[0].append(i/8)
    h[1].append((i/4)%2)
    h[2].append((i/2)%2)
    h[3].append(i%2)
h = sorted(h)

def genState():
    a = [None]*4
    a[0] = []
    a[1] = []
    a[2] = []
    a[3] = []

    for o in range(snplen):
        i = random.randint(1,14)
        a[0].append(i/8)
        a[1].append((i/4)%2)
        a[2].append((i/2)%2)
        a[3].append(i%2)
    return a

def allNeighbors(haps):
    out = []
    for i in range(len(haps)):
        for j in range(snplen):
            hapcopy = copy.deepcopy(haps)
            hapcopy[i][j] = 1 - haps[i][j]
            total = hapcopy[0][j] + hapcopy[1][j] + hapcopy[2][j] + hapcopy[3][j]
            if total < 4 and total > 0:
                out.append(hapcopy)
    return out
    
def randomPositionReads(chance,snps):
    i = random.randint(0,3)
    out = [2]*snplen
    for j in range(snplen):
        if random.randint(1,100) <= chance:
            out[j] = h[i][j]
    return out

def generateRead(length,snps):
    i = random.randint(0,3)
    readlen = min(random.randint(3*length/4,5*length/4),length)
    start = random.randint(0,snps-readlen)
    gapStart = random.randint(0,snps-1)
    gapEnd = random.randint(gapStart, min(snps-1,snps-readlen+max(0,gapStart-start)))
    out = [2]*snplen
    if gapStart > start and start+readlen > gapStart:
        for j in range(start,gapStart):
            out[j] = h[i][j]
        for j in range(gapEnd,gapEnd+readlen-gapStart+start):
            out[j] = h[i][j]
    else:
        for j in range(start,start+readlen):
            out[j] = h[i][j]
    return out

def nTest(size):
    global h
    global reads
    h = [None]*4
    
    h[0] = []
    h[1] = []
    h[2] = []
    h[3] = []

    for o in range(snplen):
        i = random.randint(1,14)
        h[0].append(i/8)
        h[1].append((i/4)%2)
        h[2].append((i/2)%2)
        h[3].append(i%2)
    h = sorted(h)
    reads = []
    for i in range(50):
        reads.append(generateRead(size,snplen))
    reads = sorted(reads)

def check(reads,guess):
    for read in reads:
        match = False
        for line in guess:
            for i in range(len(read)):
                if read[i] != 2 and read[i] != line[i]:
                    break
                elif i == len(read) - 1:
                    match = True
            if match:
                break
        if not match:
            return False
    return True

def errorCalc(reads,h):
    total = 0
    for read in reads:
        error = snplen+1
        for line in h:
            lineError = 0
            for i in range(len(read)):
                if read[i] != 2 and read[i] != line[i]:
                    lineError += 1
            error = min(lineError,error)
        total += error
    return total


def brute(reads):
    guessed = [None]*4
    guessed[0] = [0]*snplen
    guessed[1] = [0]*snplen
    guessed[2] = [0]*snplen
    guessed[3] = [1]*snplen
    while True:
        if check(reads,guessed):
            return guessed
        for i in reversed(range(snplen)):
            score = guessed[0][i]*8+guessed[1][i]+guessed[2][i]+guessed[3][i]
            if score == 14:
                guessed[0][i] = 0
                guessed[1][i] = 0
                guessed[2][i] = 0
                guessed[3][i] = 1
            else:
                score += 1
                guessed[0][i] = score/8
                guessed[1][i] = (score/4)%2
                guessed[2][i] = (score/2)%2
                guessed[3][i] = score%2
                break

def greedy(reads):
    guessed = [None]*4
    guessed[0] = [2]*snplen
    guessed[1] = [2]*snplen
    guessed[2] = [2]*snplen
    guessed[3] = [2]*snplen
    for read in reads:
        for j in range(len(guessed)):
            match = False
            for i in range(len(read)):
                if read[i] != 2 and guessed[j][i] != 2 and read[i] != guessed[j][i]:
                    break
                elif i == len(read) - 1:
                    match = True
            if match:
                for i in range(len(read)):
                    if read[i] != 2:
                        guessed[j][i] = read[i]
                break
    for i in range(len(guessed)):
        for j in range(snplen):
            if guessed[i][j] == 2:
                if guessed[0][j] != 0 and guessed[1][j] != 0 and guessed[2][j] != 0 and guessed[3][j] != 0:
                    guessed[i][j] = 0
                else:
                    guessed[i][j] = 1
    return sorted(guessed)

def hillSearch(reads):
    current = greedy(reads)
    error = errorCalc(reads,current)
    for i in range(10):
        if error == 0:
            return current
        for h in allNeighbors(current):
            h_err = errorCalc(reads,h)
            if h_err < error:
                error = h_err
                current = h
                break
    return current
    
def offOriginal(guess):
    error = snplen*4+1
    for g in list(itertools.permutations(guess)):
        part = 0
        for i in range(len(guess)):
            for j in range(snplen):
                if h[i][j] != g[i][j]:
                    part += 1
        error = min(part,error)
    return error
    
def timeTest(readLength,alg):
    totalTime = 0
    totalError = 0
    for i in range(20):
        nTest(readLength)
        start = time.clock()
        guess = alg(reads)
        totalTime += time.clock()-start
        totalError += offOriginal(guess)
    print totalError/20.0
    return totalTime/20.0
        
    