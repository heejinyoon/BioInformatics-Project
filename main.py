# D18123905 - Heejin Yoon
# Finding ORF in prokaryotic genomes
# Show the frame number
# The hucleotide and amino acid start and stop positions
# The length of the ORFs in nucleotide/ amino positions and display the DNA/AA sequence for all potential Realistic ORF.in the Sars Cov 2 virus

# 1. Open the file containing DNA nucleotide sequence (In fasta file format bp contain ~60 bp per line);
# 2. ignore descriptor line, convert each line into a single sequence (do not forget to remove all EOL markers.)
# 3. Translate the sequence ( 1st reading frame) into an amino acid sequence.
# 4. Shift one position to the right (2nd reading frame) and translate into an amino acid sequence; repeat for reading frame three.
# 5. Get the reverse compliment of the initial sequence and repeat the process (2-3)
# 6. Mark the start and the stop amino acids
# 7. Look for sequences with a start followed by a stop if there is none then there is no ORF in that reading frame.
# 8. If there is a potential ORF see if it maybe spurious using methods referred to earlier;
#    e.g.  Determine length of ORF and if less than 20 eliminate as it is a “false positive”

import sys
import tkinter
from tkinter import filedialog

def main():
    root = tkinter.Tk()
    root.wm_withdraw()  # this completely hides the root window
    # use windows explorer to input the file name
    FileName = filedialog.askopenfilename(filetypes=[('All files', '*.*')])
    root.destroy()

    # Open the file containing DNA nucleotide sequence
    DNA = Read(FileName)

    # reverse the DNA sequence
    reverse = DNA[::-1]

    # Get the reverse compliment of the initial sequence
    compliment = GetCompliment(reverse)

    # Divide reading frames - positive and negative
    # List 3 positive reading frames
    positiveKey = ["+1", "+2", "+3"]
    # List 3 negative reading frames
    negativeKey = ["-1", "-2", "-3"]
    positiveRF, negativeRF = ReadingFrame(DNA, compliment)

    # Translate the reading frames into Amino Acid Sequence
    AminoAcidList, AminoAcidSeq = DNAToAminoAcid(positiveRF, negativeRF)

    # print the 6 reading frames split into positive and negative
    print("The Complete DNA sequences obtain using concatenation is:")
    for i in range(0, len(positiveRF)):
        print("\n" + positiveKey[i] + " reading frame primary strand: \n")
        for j in range(0, len(positiveRF[i])):
            print(positiveRF[i][j])

    for i in range(0, len(negativeRF)):
        print("\n" + negativeKey[i] + " reading frame compliment strand: \n")
        for j in range(0, len(negativeRF[i])):
            print(negativeRF[i][j])

    # Mark the start and the stop amino acids, Also Find ORFs within the sequence.
    startPos, stopPos, ORF = StartStop(AminoAcidList, AminoAcidSeq)

    # Display potential ORF
    DisplayORF(AminoAcidSeq, startPos, stopPos, ORF)

def Read(FileName):
    # open the file in reading text mode using error checking
    try:
        Fp1 = open(FileName, 'r')  # Open the file containing DNA nucleotide sequence
        Fp2 = open(FileName, 'r')  # Open the file again to get descriptor line

        # read the file without description line
        Data = Fp1.readlines()[1:]

        # read the descriptor line
        descriptor = Fp2.readlines(1)

        # print output
        print("The following is the descriptor line and amino acid seq of the file:")
        print(descriptor)
        print('\n')
        print(Data)

        # convert each DNA line into a single sequence and remove all EOL markers
        FormatData = ''.join(Data)
        EndLine = FormatData.split('\n')
        DNA = ''.join(EndLine)
        print(DNA)

    except IOError:
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        stop = input()
        Fp1.close()
        sys.exit(1)

    return DNA

def ReadingFrame(DNA, compliment):
    # create list to store 3 positive reading frames - primary strand
    PositiveRF = []

    # create list to store 3 negative reading frames - complimentary strand
    NegativeRF = []

    # divide the DNA sequence into 3 lists
    for i in range(0, 3):
        PositiveRF.append([DNA[i::]])

    # divide the compliment sequence into 3 lists
    for i in range(0, 3):
        NegativeRF.append([compliment[i::]])

    return PositiveRF, NegativeRF

def GetCompliment(reverse):
    # store the compliments in the dictionary
    compDict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    reverseDna = list(reverse)

    # replace the complements with the compDict dictionary
    for i in range(len(reverseDna)):
        reverseDna[i] = compDict[reverseDna[i]]

    # make it as a string
    compString = ''.join(reverseDna)

    return compString

def DNAToAminoAcid(positiveRF, negativeRF):
    # store the codon into the array
    codon = []
    tempCodon = []

    # use a while loop to divide DNA sequence into codons
    i = 0
    while i < len(positiveRF):
        for n in range(0, len(positiveRF[i])):
            RF = positiveRF[i][n]
            for j in range(0, len(RF), 3):
                tempCodon.append(RF[j:j+3])
        codon.append(tempCodon)
        tempCodon = []
        i += 1

    # use a while loop to divide compliment DNA sequence into codons
    i = 0
    while i < len(negativeRF):
        for n in range(0, len(negativeRF[i])):
            RF = negativeRF[i][n]
            for j in range(0, len(RF), 3):
                tempCodon.append(RF[j:j + 3])
        codon.append(tempCodon)
        tempCodon = []
        i += 1

    # DNA <-> AA translation table: CodonTable
    CodonTable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    # declare an empty list an empty string
    AminoAcidList = []
    AminoAcidSeq = ''
    tempAminoAcid = []
    AA = []

    # use a nested for loop to get amino acid sequence
    for i in range(0, len(codon)):
        for j in range(0, len(codon[i])):
            if codon[i][j] in CodonTable: # translate the codon into an amino acid
                AminoAcid = CodonTable[codon[i][j]]
                AA.append(AminoAcid)
                tempAminoAcid.append(AminoAcid)
                AminoAcidSeq = ''.join(AA)
        AminoAcidList.append(tempAminoAcid)
        tempAminoAcid = []

    return AminoAcidList, AminoAcidSeq

def StartStop(AminoAcidList, AminoAcidSeq):
    # determine and display all the starts and display the positions
    # go through each amino acid and if it is start display amino acid number

    # set index to start of sequence
    index = 0
    startList = []
    stopList = []
    tempORF = []
    ORF = []

    # go through the Amino Acid sequence and find the Start and Stop point
    while index < len(AminoAcidSeq):
        if AminoAcidSeq[index] == 'M':
            startList.append(index + 1)
            tempIndex = index + 1 # increment to next AA
            while AminoAcidSeq[tempIndex] != '*' and tempIndex < len(AminoAcidSeq):
                tempIndex += 1
            storeIndex = tempIndex
            stopList.append(storeIndex + 1)
            index = tempIndex
        else:
            index += 1

    # go through the each list of Amino Acid and find the ORF
    # Look for sequences with a start followed by a stop if there is none then there is no ORF in that reading frame
    for i in range(0, len(AminoAcidList)):
        j = 0 # start point
        # when the start point found, it will keep adding to the ORF until found the stop point or it reached the end of the list
        while j < len(AminoAcidList[i]):
            if AminoAcidList[i][j] == 'M':
                tempORF.append(AminoAcidList[i][j])
                # sliding window algorithm, keep moving the position
                k = j + 1
                while AminoAcidList[i][k] != '*' and k + 1 != len(AminoAcidList[i]):
                    tempORF.append(AminoAcidList[i][k])
                    k += 1
                # Store all the ORFs that found
                ORF.append(tempORF)
                tempORF = []
                # Sliding window algorithm, when the next iteration starts, make sure it start after the stop point
                j = k
            else:
                j += 1

    return startList, stopList, ORF

def DisplayORF(AminoAcidSeq, startPos, stopPos, ORF):
    print('The Amino Acid sequence is:\n')
    print(AminoAcidSeq)

    count = 1
    # Determine length of ORF and if less than 20 eliminate as it is a “false positive”
    for i in range(0, len(startPos)):
        if len(ORF[i]) > 20:
            print('The ORF' + str(i+1) + '; The ORF start is: ' + str(startPos[i]) + '; The ORF end is: ' +
                  str(stopPos[i]) + '; The length is: ' + str(len(ORF[i])))
            ORFstr = ''.join(ORF[i])
            count += 1
            print('\n' + ORFstr)
    print('\nTotal ORFs: ' + str(count))

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

