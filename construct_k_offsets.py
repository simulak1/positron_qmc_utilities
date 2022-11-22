import argparse

def get_args():
  """Define the task arguments with the default values.
   Returns:
   experiment parameters                                                                              
  """

  args_parser = argparse.ArgumentParser()

  # Data files arguments
  args_parser.add_argument(
      '--files',
      help='twist files.',
      nargs='+',
      type=str,
      default=['1']
  )

  
  return args_parser.parse_args()


def main():
  args=get_args()
  with open("k_offsets","w") as fto:
    fto.write(str(len(args.files))+"\n")
    for j in range(len(args.files)):
      file=args.files[j]
      kpoint=1
      with open(file+"/pwfn.data") as f:
        lines=f.readlines()
        for i in range(len(lines)):
          words=lines[i].split()
          if len(words)>2 and words[0]=="k-point" and words[1]=="#":
            if(kpoint==1 and j==0):
              words=lines[i+1].split()
              fto.write(words[3]+" "+words[4]+" "+words[5]+"\n")
              break
            elif(kpoint==2 and j > 0):
              words=lines[i+1].split()
              fto.write(words[3]+" "+words[4]+" "+words[5]+"\n")
              break
            kpoint+=1  
                

if __name__ == '__main__':
    main()
