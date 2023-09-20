#!/user/bin/python3

import numpy as np
import argparse
import itertools

class smith_waterman:

    def __init__(self, seq1, seq2, match, mismatch, gaps, name):
        '''
        Parameters
        ----------
        seq1 : zstr
            The first sequence to align.
        seq2 : str
            The second sequence to align.
        match : float
            Gain in score for base match.
        mismatch : float
            Lost in score for base mismatch.
        gaps : float
            Scale of lost in score for gap.
        name : str
            Name for txt file of the output.


        Smith-Waterman class performa a local alignment by aligning two seqeunces 
        by populating self._matrix with the scores based on the values given by the user 
        and self._matrix_steps with the movements (correspondinge to the usual arrows you 
        see in every theory lesson).
        It return the bests local alignemnts, and it will be print the solution and also save them
        inside a file.txt generate, the name of the file can be change by the user
        '''
        self._first_seq = seq1.upper()  # first sequence 
        self._second_seq = seq2.upper() # second sequence
        self._value_match = match # value attribuated to mathces 
        self._value_mismatch = mismatch # value attribuated to mismatches 
        self._value_gaps = gaps # value attribuated to gaps
        self._matrix = np.zeros((len(seq1) + 1,len(seq2) + 1))  # matrix.shape [0] are the rown and matrix.shape[1] are the columns
        self._matrix_steps = [['' for col in range(len(seq2)+1)] for row in range(len(seq1)+1)] # row of the marìtrix 
        self._name = name 
        self.list_for_traceback = [''] # list used to containe the path of the traceback
        self._pointer = 0  # used in traceback to change the 
        self._dictionary_of_path = {} # containe all the path for each starter we take in consideration 

    '''Generation of the matrices, it fill the empty matrices, both self._matrix 
    and self._matrix_steps, with the scores or the paths(arrows) '''
    def matrix(self):
        first_sequence = self._first_seq
        second_sequence = self._second_seq
        match_value = self._value_match
        mismatch_value = self._value_mismatch
        gap_value = self._value_gaps
                
        for i,j in itertools.product(range(1,(self._matrix.shape[0])),range(1,(self._matrix.shape[1]))):  
                      
            if first_sequence[ i - 1] == second_sequence[ j - 1]:
                match_mismatch = (self._matrix[i - 1 , j - 1] + match_value)
            else :
                match_mismatch = (self._matrix[i - 1 , j - 1] + mismatch_value)
            
            delete = (self._matrix[i - 1, j] + gap_value) # vertically
            insert = (self._matrix[i, j - 1] + gap_value) # horizonatlly 
            best = max(match_mismatch, delete, insert,0)
           
            self._matrix[i,j] = best 
            self._matrix_steps[i][j] = self._matrix_steps[i][j] + self.genereation_matrxi_symbol([0,match_mismatch,delete,insert],best)

        # print(self._matrix_steps)
        # print(self._matrix)

    '''Creation of the matrix which containes the paths
    
        Parameters
        ----------
        list_of choices : list
                    it containe the values of the cells that are converging on the cell that is consider 
        choice : float
                    rapresent the chose number for the self._matrix which need to take into consideration

        So based on the values the ones that mathc the best value is possible to undertstend the path which are coming 
        
        '''
    def genereation_matrxi_symbol(self,list_of_choices,choice):
        momentary = []
        result =''

        for i,x in enumerate(list_of_choices): 
            if choice == x:                  # i is the index and x is the value of the movements 
                momentary.append(i)

        if len(momentary) == 3:  # all three ( L = left, D = diagonal, U = up)
            result = result + 'LDU'
        elif len(momentary) == 1:
            if momentary[0] == 0: # N means stops
                result = result + '' # value zero so stop 
            elif momentary[0] == 1:
                result = result + 'D' # match or mismatch 
            elif momentary[0] == 2:
                result = result + 'U' #upward
            elif momentary[0] == 3:
                result = result + 'L' # leftward
        elif len(momentary) == 2: # based on the sums 
            if (momentary[0] + momentary[1]) == 3:
                result = result + 'DU'  #match and insert, match - upwatd 
            elif (momentary[0] + momentary[1]) == 4:
                result = result + 'DL' # match and deletion, match - leftward
            elif (momentary[0] + momentary[1]) == 5:
                result = result + 'LU' # delete and insert, upward - leftward 
        return result
    
    '''Section to insert each coordinates in the backtracking section and give the results of the score function
        it's a sort of 'linker' function of the next fucntions '''
    def alignment(self):
        scoring_matrix = self._matrix
        coord_higher_values = self.best_scores(scoring_matrix)
        # print(coord_higher_values)

        for i in coord_higher_values:
            # print(i)
            self.list_for_traceback = ['']
            self._pointer = 0 
            self.moves_and_scores(i,i)
        # print(self._dictionary_of_path)
        
        return self.scores()
    
    '''Function where retun the coordinates of the highest scores
        Parameters
        ----------
        
        scoring_matrix = it's the matrix that containe all the scores 
         
     '''
    def best_scores(self,scoring_matrix): # change name into biggest-scores
        
        coordinates = []
        highest_value = np.amax(scoring_matrix)     # find the maximum values
        coordinates_raw = np.where(scoring_matrix == highest_value)
        
        for i in range(len(coordinates_raw[0])):
            coordinates.append((coordinates_raw[0][i],coordinates_raw[1][i]))
        return coordinates

    '''Section where happen the backtracking by saving each path based on the matrix whihc containe the steps
       So by starting from the best score/scores it make a backtrack base on the self._matrix_step and saved into dictionary_of_path
       
       Parameters
       ----------
        coordinates_start : tuple
                        it contain the coordinates where the cell of interest is, it change time get recursed the fucntion
        coordinates_keys : tuple
                        containe the coordinates of the start position of the best score, it remain constant and they will be the keys for the dictionary_of_path
       
       '''
    def moves_and_scores(self,coordinates_start,coordinates_keys):
        i,j = coordinates_start[0],coordinates_start[1]
        cell = self._matrix_steps[i][j]
        # print(i,j, cell,len(cell))
        if len(cell) == 1:
            if cell == 'D':
                self.list_for_traceback[self._pointer] += 'D'
                # print(self.list_for_traceback)
                self.moves_and_scores((i-1,j-1),coordinates_keys)
            elif cell=="L":
                self.list_for_traceback[self._pointer]+="L"
                # print(self.list_for_traceback)
                self.moves_and_scores((i, j-1),coordinates_keys)
            elif cell=="U":   
                self.list_for_traceback[self._pointer]+="U"
                # print(self.list_for_traceback)
                self.moves_and_scores((i-1, j),coordinates_keys)
        elif cell == '':
            # print("cella vuota, fine")
            if coordinates_keys not in self._dictionary_of_path:
                # print(coordinates_keys)
                self._dictionary_of_path[coordinates_keys]=[self.list_for_traceback[self._pointer]]
            else:
                self._dictionary_of_path[coordinates_keys].append(self.list_for_traceback[self._pointer])
            self.list_for_traceback.remove(self.list_for_traceback[self._pointer])
            self._pointer=len(self.list_for_traceback)-1
            # print(self.list_for_traceback)
            # print(self._pointer)
        
        elif len(cell) > 1:
            # print(cell, "Len cell:", len(box))
            temporary = []
            temporary.append((self.list_for_traceback[self._pointer]))
            # print("temporary",temporary, self.list_for_traceback[self._pointer])
            self.list_for_traceback.extend(temporary * (len(cell)-1) )
            # print(self.list_for_traceback, "lunghezza listina", len(self.list_for_traceback))
            
            for t in cell:
                if t =="D":                  
                    self.list_for_traceback[self._pointer]+="D"
                    # print(self.list_for_traceback)
                    self.moves_and_scores((i-1, j-1),coordinates_keys)
                elif t =="L":
                    self.list_for_traceback[self._pointer]+="L"
                    # print(self.list_for_traceback)
                    self.moves_and_scores((i, j-1),coordinates_keys)
                elif t=="U":
                    self.list_for_traceback[self._pointer]+="U"
                    # print(self.list_for_traceback)
                    self.moves_and_scores((i-1, j),coordinates_keys)
            # print(self.list_for_traceback)
            # print(self._pointer)
        # print("duct",self._dictionary_of_path)
       
    '''Section for the part of scores based on the path we have inside the dictionary.
       So base on the self.dicitonary_of_path which contain the path of each alignemnt with the best score
       the function take the path information and gnerate the scores,the part of seqeunces that are aligned 
       and the symbols that is possible to see during the printing 
    '''
    def scores(self):
        dictionary_with_alignment = {}
        for i in self._dictionary_of_path.keys():
            # print(self._dictionary_of_path)
            for j in self._dictionary_of_path[i]:
                seq_1 = ''
                seq_2 = ''
                symbols = ''
                n_match = 0 
                n_mismatch = 0
                n_gap = 0
                lenght = 0
                x,y = i[0]-1,i[1]-1
                for k in j:
                    # print("symbol:", k, 'seqeunces:',self._first_seq[x],self._second_seq[y])
                    if k == 'D':
                        if self._first_seq[x] == self._second_seq[y]:
                            seq_1 += self._first_seq[x]
                            seq_2 += self._second_seq[y]
                            symbols += '*'
                            n_match += 1 
                            lenght += 1
                        else:
                            seq_1 += self._first_seq[x]
                            seq_2 += self._second_seq[y]
                            symbols += '|'
                            n_mismatch += 1 
                            lenght += 1
                        x,y = x-1, y-1 
                    elif k == 'U':
                        seq_1 += self._first_seq[x]
                        seq_2 += '_'
                        symbols += ' '
                        n_gap += 1 
                        lenght += 1
                        x = x-1
                    else:
                        seq_2 += self._second_seq[y]
                        seq_1 += '_'
                        symbols += ' '
                        n_gap += 1 
                        lenght += 1
                        y = y -1
                listina_scores = [seq_1[::-1],seq_2[::-1],symbols[::-1],n_match,n_mismatch,n_gap,lenght]
                if self._matrix[i] not in dictionary_with_alignment:
                    dictionary_with_alignment[self._matrix[i]] = []
                    dictionary_with_alignment[self._matrix[i]].append(listina_scores)
                else:
                    dictionary_with_alignment[self._matrix[i]].append(listina_scores)
        # print(dictionary_with_alignment)

        return dictionary_with_alignment

    ''' section part where happen the print part of the alignments and generation and wiriting of them in the file txt  '''
    def printer(self): 
        
        self.matrix()
        dictionary_to_print = self.alignment()
        # print(dictionary_to_print)

        content = ''
        m = 0

        for i in dictionary_to_print:
            selcted = dictionary_to_print[i]
            for k in selcted:
                m += 1
                titles = '\nOutput {}:'.format(m)
                final_setting = titles + '''\nThe score of the alignemnt:{}\tlength of the alignment:{}\tnumber of match :{}\tnumber of mismatch:{}\tnumber of gaps:{}\n{}\n{}\n{}
                '''.format(i,k[6],k[3],k[4],k[5],k[0],k[2],k[1])
                print(final_setting)
                content += final_setting

        f = open(self._name,'w') # part of creating the txt file of the alignment 
        f.write(content)
        f.close()
                
if __name__ == '__main__':
            
    parser = argparse.ArgumentParser(description='''Smith-Waterman algorithm is for optimal local sequence alignment. Use as input sequences by typing them manually.
    \nExample : python3 smith-waterman.py ATCGGCGATA ATTATACGATA -g -2 -m 3 -p -1 -o alignment_output''', epilog='The output will be written inside a file.txt, which will be generated by running the algorithm in the same folder of the file code.',formatter_class= argparse.RawTextHelpFormatter)
    # use -h to get the helper, not needed to explicit it, there is already in default
    parser.add_argument('seq1', type=str, help='First sequence, which correspond of the matrix rows')
    parser.add_argument('seq2', type=str, help='Second sequence, which correspond of the matrix columns')
    parser.add_argument('-g','--gap_penalty', type=float, default= -1, help='Scoring value of a gap. Default is -1')
    parser.add_argument('-m','--match', type=float, default=2, help='Scoring value of a match. Default is 2')
    parser.add_argument('-p','--mismatch', type=float, default= -1, help='Scoring value of a mismatch. Default is -1')
    parser.add_argument('-o','--output_name', type=str, default='SW_output',help='Name of the output file that will be created (default: SW_output)')
    args = parser.parse_args()

    name_out= args.output_name + '.txt'
    '''To print the name of the output'''
    print('The name of the output is: {}'.format(name_out))

    sw = smith_waterman(args.seq1, args.seq2,gaps=args.gap_penalty, match=args.match, mismatch=args.mismatch, name=name_out)
    sw.printer()
    
