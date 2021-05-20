file = open('P:/CodingProject/BioInformatics/BioInformatics/Alignment Sequences/274326_GlobAlighn_Sequences.txt','r')

file_contents = file.readlines()

while file_contents:

    s1=file_contents[0][6:23]
    print("This is Sqeuence 1 : ", s1)
    leng1 = len(s1)
    
    s2=file_contents[1][6:23]
    print("This is Sqeuence 2 : ", s2)
    leng2 = len(s2)

    print("\n")
    file_contents=file.readline()

matrix  = []

gap = 0
match = 1
mismatch = 0

def matrixTable(score,s1,s2):

    s1 = "   " + s1;
    s2 = " " + s2;

    for i in range(0, len(s1)):
        print(" ",s1[i], end = "\t")

    for i in range(0, len(score)):
        print("\n", s2[i], end = "\t")
        print(" ", end = "\t")
        for j in range(0, len(score[i])):
            print("[", score[i][j], "]", end = "\t")
        print(" ")

def zeros(rows, cols):

    retval=[]

    for i in range(rows):
        retval .append([])
        for j in range(cols):
            retval [-1].append(0)
    return retval

def checkingMatchOrMismatch(seq1, seq2):
    if seq1 == seq2:
        return match
    else:
        return mismatch

def GlobalAlignment(s1, s2):

    score = zeros(leng1+1, leng2+1)

    for i in range(0, leng1+1):
        score[i][0] = gap * i

    for j in range(0, leng2+1):
        score[0][j] = gap * j

    for i in range(1, leng1+1):
        for j in range(1, leng2+1):
            diagonal = score[i-1][j-1] + checkingMatchOrMismatch(s1[j-1], s2[i-1])
            top = score[i-1][j] + gap
            left = score[i][j-1] + gap

            score[i][j] = max(diagonal, top, left)
    return score

matrixTable(GlobalAlignment(s1, s2),s1,s2)

def traceBack(s1, s2):

    score = GlobalAlignment(s1, s2)
    output1 = ""
    output2 = ""
    outputSymbol = ""

    count = 0
    gap = 0
    match = 1
    mismatch = 0

    i = leng2
    j = leng1

    
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_top = score[i][j-1]
        score_left = score[i-1][j]
        
        if score_current == score_diagonal + checkingMatchOrMismatch(s1[j-1], s2[i-1]):
            if(s1[j-1] == s2[i-1]):
                count += match
                output1 += s1[j-1]
                outputSymbol += "|"
                output2 += s2[i-1]
            else:
                count+=mismatch
                output1 += s1[j-1]
                outputSymbol += " "
                output2 += s2[i-1]
            i -= 1
            j -= 1
        elif score_current == score_top + gap:
            count += gap
            output1 += s1[j-1]
            outputSymbol += " "
            output2 += '-'
            j -= 1
        elif score_current == score_left + gap:
            count += gap
            output1 += '-'
            outputSymbol += " "
            output2 += s2[i-1]
            i -= 1

    while j > 0:
        output1 += s1[j-1]
        outputSymbol += " "
        output2 += '-'
        j -= 1
    while i > 0:
        output1 += '-'
        outputSymbol += " "
        output2 += s2[i-1]
        i -= 1
    
    output1 = output1[::-1]
    outputSymbol = outputSymbol[::-1]
    output2 = output2[::-1]
    
    return(output1, outputSymbol, output2, count)

output1, outputSymbol, output2, count = traceBack(s1, s2)
print("\nTrace back Result: \n", output1 + "\n " + outputSymbol + "\n " + output2, "\nOptimal Trace back Result: ", count, "\n")

file.close()