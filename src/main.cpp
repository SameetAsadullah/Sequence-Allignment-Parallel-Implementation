#include <mpi.h>
#include <omp.h>
#include <vector>
#include <fstream>
#include <string.h>
#include <iostream>
using namespace std;

// basic global variables
int w1, w2, w3, w4, c1_size = 9, c2_size = 11;
string conservative_grp[] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
string semi_conservative_grp[] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};


// functions for mutation
bool checkSubstitution(char first, char second) 
{
    for (int i = 0; i < 9; ++i) 
    {
        if (conservative_grp[i].find(first) != std::string::npos && conservative_grp[i].find(second) != std::string::npos) 
        {
            return false;
        }
    }
    
    return true;
}

vector<string> getMutants(string seq) 
{
    vector<string> mutants;
    
    #pragma omp parallel for schedule(static, 1), collapse(2)
    for (int i = 0; i < seq.length(); ++i) 
    {
        for (int j = 0; j < 26; ++j) 
        {
            if (seq[i] != char(j + 65) && checkSubstitution(seq[i], char(j + 65))) 
            {
                string seq_copy = seq;
                seq_copy[i] = char(j + 65);
                
                #pragma omp critical
                {
                	mutants.push_back(seq_copy);
                }
            }
        }
    }
    return mutants;
}


// function to calculate alignment score
int stringComparison(string seq1, string seq2)
{
	bool flag = 0;
	int no_of_stars = 0, no_of_colons = 0, no_of_dots = 0, no_of_spaces = 0;
	
	for(int i = 0; seq2[i] != '\0'; i++)
	{
		if (seq1[i] == seq2[i])
		{
			no_of_stars++;
			continue;
		}
		
		for (int j = 0; j < c1_size; j++)
		{
			if (conservative_grp[j].find(seq1[i]) != string::npos && conservative_grp[j].find(seq2[i]) != string::npos)
			{
				no_of_colons++;
				flag = 1;
				break;
			}
		}
		
		if (flag == 1)
		{
			flag = 0;
			continue;
		}
		
		for (int j = 0; j < c2_size; j++)
		{
			if (semi_conservative_grp[j].find(seq1[i]) != string::npos && semi_conservative_grp[j].find(seq2[i]) != string::npos)
			{
				no_of_dots++;
				flag = 1;
				break;
			}
		}
		
		if (flag == 1)
		{
			flag = 0;
			continue;
		}
			
		no_of_spaces++;
	}
	
	return w1*no_of_stars - w2*no_of_colons - w3*no_of_dots - w4*no_of_spaces;
}

int main()
{
	// MPI variables and initialization
	int size, rank, root = 0;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// to store information extracted from input file
	bool minmax;
	string seq1, seq2;
	
	// obtaining values for w1, w2, w3 and w4 globally
	string line;
	fstream file("input.txt", ios::in);
	
	file >> line;
	w1 = stoi(line);
	file >> line;
	w2 = stoi(line);
	file >> line;
	w3 = stoi(line);
	file >> line;
	w4 = stoi(line);
	
	getline(file, seq1);
	getline(file, seq1);
	getline(file, seq2);
	getline(file, line);
	
	if (line == "minimum")
		minmax = 0;
	else
		minmax = 1;
	
	file.close();
	
	vector <string> vec_str = getMutants(seq2);
	
	if (rank == root)
	{	
		// MPI task division
		if (seq1.length() > seq2.length())
		{
			line = seq1;
			int slave = 1, score, store, offset = 0, offset_store;
			bool check = 0, flag = 0;
			char to_recv[5000];
			string str_store;
			
			for(int i = 1; line.size() >= seq2.size(); i++)
			{
				MPI_Send(line.c_str(), line.length(), MPI_CHAR, slave, 1, MPI_COMM_WORLD);
				
				line = seq1.substr(i, seq1.size() - 1);
				
				if (check == 1)
				{
					MPI_Recv(&score, 1, MPI_INT, slave, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&to_recv, 5000, MPI_CHAR, slave, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					
					if (flag == 1)
					{
						store = score;
						str_store = to_recv;
						offset_store = offset;
					}
					else if (minmax == 0 && score < store)
					{
						store = score;
						str_store = to_recv;
						offset_store = offset;
					}
					else if (minmax == 1 && score > store)
					{
						store = score;
						str_store = to_recv;
						offset_store = offset;
					}
					
					flag = 0;
					offset++;
				}
				
				slave++;
				if (slave >= size)
				{
					if (check == 0)
						flag = 1;
					
					slave = 1;
					check = 1;
				}
			}
			
			if ((size - 1) > (seq1.size() - seq2.size() + 1))
			{
				slave = 1;
				size = seq1.size() - seq2.size() + 2;
			}
			
			for (int i = 1; i < size; i++)
			{
				MPI_Recv(&score, 1, MPI_INT, slave, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&to_recv, 5000, MPI_CHAR, slave, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				if (minmax == 0 && score < store)
				{
					store = score;
					str_store = to_recv;
					offset_store = offset;
				}
				else if (minmax == 1 && score > store)
				{
					store = score;
					str_store = to_recv;
					offset_store = offset;
				}
				
				slave++;
				offset++;
				if (slave >= size)
				{
					slave = 1;
					break;
				}
			}
			
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			
			line = "abort";
			for(int proc = 1; proc < size; proc++)
			{
				MPI_Send(line.c_str(), line.length(), MPI_CHAR, proc, 1, MPI_COMM_WORLD);
			}
			
			fstream ofile("output.txt", ios::out);
			ofile << str_store << endl;
			ofile << offset_store << " " << store;
			ofile.close();
		}
		
		else if (seq1.length() == seq2.length())
		{
			line = "abort";
			for(int proc = 1; proc < size; proc++)
			{
				MPI_Send(line.c_str(), line.length(), MPI_CHAR, proc, 1, MPI_COMM_WORLD);
			}

			int score = stringComparison(seq1, seq2);
			int store = score;
			string to_send = seq2;
			
			for (int i = 0; i < vec_str.size(); i++)
			{
				score = stringComparison(seq1, vec_str[i]);
				
				if (minmax == 0 && score < store)
				{
					store = score;
					to_send = vec_str[i];
				}
				else if (minmax == 1 && score > store)
				{
					store = score;
					to_send = vec_str[i];
				}
			}
			
			fstream ofile("output.txt", ios::out);
			ofile << to_send << endl;
			ofile << "0 " << store;
			ofile.close();
		}
		
		else
		{
			line = "abort";
			for(int proc = 1; proc < size; proc++)
			{
				MPI_Send(line.c_str(), line.length(), MPI_CHAR, proc, 1, MPI_COMM_WORLD);
			}
			
			cout << "ERROR: Sequence 1 is shorter than Sequence 2!" << endl;
			cout << "Program terminating..." << endl;
		}
	}
	
	else
	{
		char recv1[10000];
		
		while(true)
		{
			MPI_Recv(recv1, 10000, MPI_CHAR, root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			if (strcmp(recv1, "abort") == 0)
			{
				break;
			}
			
			string str_recv = recv1;
			int score = stringComparison(str_recv, seq2);
			int store = score;
			string to_send = seq2;
			
			#pragma omp parallel for schedule(static, 1)
			for (int i = 0; i < vec_str.size(); i++)
			{
				int temp = stringComparison(str_recv, vec_str[i]);

				#pragma omp critical
				{
					score = temp;
					if (minmax == 0 && score < store)
					{
						store = score;
						to_send = vec_str[i];
					}
					else if (minmax == 1 && score > store)
					{
						store = score;
						to_send = vec_str[i];
					}
				}
			}
			
			MPI_Send(&store, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
			MPI_Send(to_send.c_str(), to_send.length(), MPI_CHAR, 0, 3, MPI_COMM_WORLD);
			
			memset(recv1, 0, sizeof recv1);
		}
	}
	
	MPI_Finalize();
}
