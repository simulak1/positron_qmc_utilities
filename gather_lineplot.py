import numpy as np
import argparse

'''
This script gathers lineplot.dat-files constructed by post-processing
expval.data-files that are the expectation value files of CASINO. The script
gathers lineplot.dat-files having two columns of length N into a numpy matrix
of dimensions (N,number_of_files). This numpy matrix is then saved. If asked by
the user, also the distance vector is saved. 
'''

def get_args():
  """Define the task arguments with the default values.
   Returns: 
   experiment parameters
  """

  args_parser = argparse.ArgumentParser()

  args_parser.add_argument(
      '--files',
      help='Folder names from which to gather lineplots..',
      nargs='+',
      type=str,
      default=['1']
  )

  args_parser.add_argument(
      '--write_r',
      help='Write a file for the distances?',
      type=int,
      default=0
  )

  args_parser.add_argument(
      '--gfile',
      help='The name of the target file',
      type=str,
      default='g'
  )
  

  return args_parser.parse_args()

def main():


    args=get_args()

    i=0
    for file in args.files:
        i+=1
        with open(file+'/lineplot.dat') as f:
            lines=f.readlines()
            N=len(lines)
            if(i==1):
                g20=np.zeros((N,len(args.files)))
                r=np.zeros((N,))
            for j in range(N):
                line=lines[j].split()
                r[j]=float(line[0])
                g20[j,i-1]=float(line[1])

    print(g20.shape)
    np.save(args.gfile,g20)
    if args.write_r>0:
        np.save('r',r)

if __name__=='__main__':
    main()
