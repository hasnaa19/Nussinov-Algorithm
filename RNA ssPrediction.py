def bifurcation(matrix, i, j):                              #This function takes the matrix, i and j values
    max_val = -1
    k_value = -1
    for k in range(i+1, j):                                 #iterate for each K value found in range i < K < j and calculate "matrix[i][k] + matrix[k+1][j]"
        if matrix[i][k] + matrix[k+1][j] > max_val:         #if the current k produces value greater than the last max_val
            max_val = matrix[i][k] + matrix[k+1][j]         #update the max value
            k_value = k                                     #update the k that produces the maximum value
    return max_val, k_value                                 #return the max value and k that produces it

############################

def are_paired(base1, base2):                                                      #This function take two bases
    if (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):         #check if they can pair with each other
        return 1                                                                   #if they can it returns 1 which is the pairing score
    elif (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A'):
        return 1
    else:
        return 0                                                                   #else it return 0 which is the unpairing score

#############################

def traceback(matrix, result, i, j, seq):                                         #this function takes the matrix,the list that will carry 2D structure of the sequence, the row/col indices at which traceback will start and the sequence
                                                                                   #in the first iteration, traceback starts at the top right corner therefore i = 0 & j = length-1
    while i < j:                                                                   #stop when??   i>j , which means that we've reached the main diagonal where i and j are equal
        paired = are_paired(seq[i], seq[j])                                        #calculate the value coming from the diagonal, left, down and bifurcation and check the current cell equals to which one of them
        left = matrix[i][j-1]
        down = matrix[i+1][j]
        diagonal = matrix[i+1][j-1] + paired
        bifurcation_val, k = bifurcation(matrix, i, j)

        if matrix[i][j] == diagonal:                                                #if we are coming from the diagonal we have two cases, paired or unpaired
            if paired:                                                              #only when we are coming from the diagonal AND two bases are paired change the result list that carries 2D structure
                #print("from diagonal paired")
                result[i] = "("                                                     #add an open bracket at the smaller between i & j and add the closed bracket at the other index, do i need to check who's the smaller?? no as above main diagonal i is always smaller that j
                result[j] = ")"
            #else:
                #print("from diagonal not paired")
            i+=1                                                                    #whether they are paired or not we will increment i, decrement j in order to move diagonally to the next cell
            j-=1
            #print(result)
        elif matrix[i][j] == down:                                                  #if the current cell's value is coming from the left, move to the cell below the current cell by incrementing rows(i) while column(j) remains the same
            #print("from down")
            i += 1
        elif matrix[i][j] == left:                                                  #if the current cell's value is coming from the left, move to the cell at the left of the current cell by incrementing column(j) while row(i) remains the same
            #print("from left")
            j-=1
        elif matrix[i][j] == bifurcation_val:                                       #if the current cell's value is coming from bifurcation at k
            traceback(matrix, result, i, k, seq)                                    #jump to matrix[i][k] and continue tracing
            traceback(matrix, result, k+1, j, seq)                                  #then jump to matrix[k+1][j] and continue tracing
            i=0                                                                     #set i & j to zero after coming back from the recursive call in order to break from the loop
            j=0

###################################

def RNA_2Dstructure_prediction(sequence):                         #This function takes the sequence that we want to predict it's 2D structure, then returns the matrix and the 2D structure
    length = len(sequence)
    result = ["." for i in range(length)]                         #initializing the list that will carry 2D structure with dots as a start

    #matrix initialization
    Matrix =[[0 for i in range(length)]for j in range(length)]    #initializing the matrix with all zeros as we gonna fill the matrix starting after the main diagonal

    #filling matrix
    for x in range(1, length):                                 #x represents the diagonal we are currently filling, since we will start filling our matrix from the diagonal right after the main one, number of iteration will be length of sequence -1
        i = 0                                                  #i represent rows while j represents columns, regardless which diagonal we are filling, we will always start from the first row
        j = x                                                  #but this is not the case with columns as in the first iteration we skipped the first column as we are not gonna fill the main diagonal ( j =1 ), and in every iteration afterward we skip a column from the start in order to move to the next diagonal
        while i < (length-x) and j < length:                   #move diagonally on each cell and stop when i < length -x , why? as we've said we will always start from the first row but we will terminate one row earlier in each iteration as a result of moving diagonally

            score = are_paired(sequence[i], sequence[j])       #check if the 2 bases corresponding to the current cell could pair with each other, if yes score =1 else score = 0
            bifurcation_val, k = bifurcation(Matrix, i, j)     #calculate the score of bifurcation if it occurred at the current cell
            Matrix[i][j] = max(                                #the value of the current cell is the max of these four choices
                Matrix[i][j-1],
                Matrix[i+1][j],
                Matrix[i+1][j-1]+score,
                bifurcation_val
            )
            # if Matrix[i][j] == bifurcation_val and Matrix[i][j] !=  Matrix[i][j-1]and Matrix[i][j] != Matrix[i+1][j] and Matrix[i][j] != Matrix[i+1][j-1]+score:
            #     print("bifurcation", bifurcation_val)
            #     print(i,j)
            i+=1                                              #increment i and j to move diagonally to the next cell
            j+=1
    traceback(Matrix, result,0 ,length-1, sequence)          #now the matrix is ready for tracing, starting the trace back from the upper right corner of the matrix which is cell [0][length-1]
    result = ''.join(x for x in result)                       #join the elements of the list in order to appear as a string
    return Matrix, result                                     #return the result of the traceback (which is one possible 2D structure for that sequence) and the matrix

#####################################

def run_algorithm(seq):                                       #this is function takes the sequence as an input and send it to RNA_2Dstructure_prediction then print the result
    matrix, result = RNA_2Dstructure_prediction(seq)
    print("MATRIX :")
    print("--------")
    for i in range(len(seq)):
        print(matrix[i])
    print("============================================")
    print("SEQUENCE:  ", seq)
    print("STRUCTURE: ", result)


#------------------------------------------MAIN-------------------------------------------

run_algorithm("GGGAAAUCC")
print("")
run_algorithm("CGGACCCAGACUUUC")

