#The following Python script is a general-purpose file list renamer.
#Frequently in this project, we must create lists of files and then append some
#prefix to each filename in the list. This script does it for us.
#
#Input:
#  - The name of the filename list file (a .txt document which lists filenames).
#  - The desired name of the output file.
#  - The prefix to be added
import sys

print("\n--- File List Renamer ---\n")

if (sys.version_info[0] < 3):                            #Aborts the script if Python version is less than 3.
  print("Error: Requires python 3. To run, type:\n") 
  print("python3 "+sys.argv[0]+"\n")
  quit()

input_filename = input("What is the name of the input file? ./") 
input_file = open(input_filename, "r")

output_filename = input("What will the output file be called? ./")
output_file = open(output_filename, "w")

prefix = input("What is the prefix to add? ")

input_lines = input_file.readlines()
output_lines = [prefix+line for line in input_lines]

input_file.close()

for line in output_lines:
  output_file.write(line)

output_file.close()
