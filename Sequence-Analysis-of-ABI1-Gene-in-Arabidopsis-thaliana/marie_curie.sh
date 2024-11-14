# Create new folder named marie_curie
mkdir marie_curie

# Create another new folder named biocomputing and redirect to it
mkdir biocomputing && cd biocomputing

# Download the working files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna 
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk 
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# Move the file format fna to the folder marie_curie
mv wildtype.fna ../marie_curie

# Remove the duplication 
rm wildtype.gbk.1

# Heading to marie_curi folder again 
cd ..
cd marie_curie

# Check if there is mutant pattern and if so, will copy Mutant
if grep -q "tatatata" "wildtype.fna"; then echo "Mutant"
fi 

# Extract the lines containing the pattern and saving them in a new file called mutant.fna
grep -n "tatatata" "wildtype.fna" > mutant.fna

# Download fasta format file to the gene ABI1 directly from NCBI
# Install EDirect software
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
# Set the path for the current terminal session 
export PATH=${HOME}/edirect:${PATH}
# Download fasta format file to the gene ABI1 using its id
efetch -db nucleotide -id X77116.1 -format fasta > ABI1.fasta

# Count the file lines except the header
grep -v '^>' ABI1.fasta | wc -l
# Count 'A' in the sequence
A_count=$(grep -v '^>' ABI1.fasta | grep -o 'A' | wc -l)
echo "A_count: $A_count"
# Count 'G' in the sequence
G_count=$(grep -v '^>' ABI1.fasta| grep -o 'G' | wc -l)
echo "G_count: $G_count"
# Count 'C' in the sequence
C_count=$(grep -v '^>' ABI1.fasta| grep -o 'C' | wc -l)
echo "C_count: $C_count"
# Count 'T' in the sequence
T_count=$(grep -v '^>' ABI1.fasta| grep -o 'T' | wc -l)
echo "T_count: $T_count"

# Calculate GC% using emboss 
file="ABI1.fasta"
total_nucleotides=$(grep -v '^>' "$file" | tr -d '\n' | wc -m)
GC_percentage=$((($G_count + $C_count) * 100 / $total_nucleotides))
echo "GC percentage: $GC_percentage %"
# Another solution to get the GC% 
source activate emboss
cat ABI1.fasta | infoseq -auto -only -pgc stdin
conda deactivate 

# Create a new fasta formate file 
touch marie_curie.fasta
# Adding and concatenate the number of 'A', 'G', 'T', and 'C'.
echo "Number of A: $A_count" >> marie_curie.fasta
echo "Number of G: $G_count" >> marie_curie.fasta
echo "Number of T: $T_count" >> marie_curie.fasta  
echo "Number of C: $C_count" >> marie_curie.fasta

# Clear the terimal and print all the commenda used 
clear
history 

# List the files in the two folder 
ls -p marie_curie
ls -p biocomputing
