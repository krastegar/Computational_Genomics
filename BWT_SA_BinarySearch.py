#-------------------------------------------------------------------------------------------------------
# Purpose: To perform a Burrow wheeler transform and suffix array on FASTA file. This is to demonstrate
#          knowledge of algorithms that use these techniques for searching for genese. Suffix Array is then 
#          used for binary search algorithim to look for local gene inside FASTA file. 
#-------------------------------------------------------------------------------------------------------
# Algorithm: 
#   Step 1: Open document as a read file and create an array that contains the entire string of nucleotides 
#           in fasta file. 
#   Step 2: I have to create an array of strings that have the last character being rotated to the beginning 
#           of the string for the entire fasta file. At the same time create an array for indicies as well 
#   Step 3: create a list of lists from the index values and strings from bwt
#   Step 4: Append each index of the list of list to a new array and then sort it (this my completed suffix array)
#   Step 5: Create the binary search algorithm that initially cuts suffix array in half and then looks at both high
#           and low ends of the data structure. 
#   Step 6: Call Binary search Algorithm recursively. If gene could not be found return error message
#-------------------------------------------------------------------------------------------------------

import numpy as np

with open("referencegenome.fasta", "r") as f:
    with open("test.txt", "w") as t: 
        array = [] # have to put array outside of the loop or else I will get a million arrays
        for line in f: # have to read the fasta in line segments otherwise i get textwrapper error
            if line.startswith(">"): # want to skip the lines that have the metadata
                pass
            else: 
                string = "".join(line.split()) # I am removing all the spaces and joining every line in fasta file together
                for nt in string: # now I can loop through this and append each nt into an array so that I can make long_string
                    array.append(nt) # appending each nt to a list without the title
        long_string ="".join(array) #joining the elements of my list/ array and turning them into a string
                                    #"" arguement joins them w/o anyspaces or symbols
        dct = []
        ind = []
        for i in range(len(long_string)):
            val = long_string[-1] + long_string[:-1] # moving first string 
            long_string = val # this updates genome so that we may continue the process of moving first letter to the front and joining the rest of the body
            dct.append(long_string) # going to build a list of a list 
            ind.append(i)
        index=np.lexsort((dct,ind)) # creating my index's for my list of list
        sa = []
        for i in range(len(index)): # have to loop throught the length of my index's to create a list of list with same length
            sa.append([dct[i],ind[i]])# using this syntax to build my list of list
        sort_bwt = sorted(sa) #sorting my bwt
        sa1=[] # 1st column of my list of list (traditionally called suffix array from teachers)
        bwt1=[] #last column of each list within my list
        for jj in sort_bwt:
            sa1.append(jj[0][0]) # index[0] is one long string index[0][0] first element of long_string
            bwt1.append(jj[-2][-2]) # I have to use [-2] becuase [-1] is the index (int) number and can't be appended and we do not want to anyways
        t.write("Suffix Array: \n" + str(sa1)+ "\n"+"BWT: \n"+str(bwt1) + "\n")
        
        def Bisect(search, ref, low, high):# low starts at 0 and high is usually len(ref)-1
            if high >= low:
                mid = (low+high)//2
                x = ref[mid][0][0:len(search)] # what this is saying is that we are looking at the first element of the first list(from 0 position all the way to len(search)) in the entire list
                if x==search: #says if our index of suffix array produces the same string as our search genome return mid (which is th midpoint of our dictionary)
                    return mid
                elif x>search: # if we are too high on our dictionary page we need to look lower so that is why our low = mid-1 and the rest stays the same
                    return Bisect(search, ref, low, mid-1)
                else: #we dont have to right ref[mid] < search because that is our only condition left that we havent tested
                # We use low=mid+1 because becuase if we are less than our search that means we have to a higher index to look for search
                    return Bisect(search, ref, mid+1, high)
            else: # since we have exhausted all options we provide an output message explaining that we could not find it
                return (-1)
        #Calling my function
        call = Bisect("GTTCCGAGAGCTGAA", sort_bwt, 0, len(sort_bwt)-1)
        # outputing results from my function call
        if call != -1:
            t.write("Bisection Algorithim: \n"+"Gene is at index:"+str(call))
        else:
            t.write("Gene is not present in reference genome")
        print("All information is on a .txt document called test.txt in your directory!")

        
