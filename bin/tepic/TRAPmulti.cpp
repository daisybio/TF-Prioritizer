//Helge Roider
//Max Planck Institute for Molecular Genetics - Berlin
//DATE: 16/06/2008

//Adapted by Florian Schmidt and Marcel H. Schulz
//Cluster of Excellence on Multimodal Computing and Interaction and Max Planck Insitute for Informatics - Saarbruecken
//DATE: 17/11/2016

//Modified by Markus Hoffmann
//DATE 07/01/2022

#include<stdio.h>
#include <stdlib.h>
#include<cmath>
#include <cfloat>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<string.h>
#include<iomanip>
#include<omp.h>
#include <vector>

using namespace std;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;


int main(int argc, char *argv[]){

    if(argc < 5){
        cerr << "\nINPUT:\n\t1. ENERGY MATRIX (from PSCM_to_PSEM)\n\t2. FASTA FILE\n\t3. NUMBER OF CORES TO BE USED\n\t4. OUTPUT SEQUENCES FILE DIRECTORY\n\t5. optional [N, Psi <default>]\n\nOUTPUT\n\tN  or  Psi (= N/seqlength)\n\n";
        exit(1);
    }
    unsigned int num_workers=atoi(argv[3]);
    unsigned int outputtype = 0;
    if(argc == 6){
        if(strcmp(argv[5],"N") == 0){
            outputtype = 1;
        }
    }

    //see argvs
    std::string argument0 = argv[0];
    std::string argument1 = argv[1];
    std::string argument2 = argv[2];
    std::string argument3 = argv[3];
    std::string argument4 = argv[4];


    //----------------------------------------------------------------
    //Initialise VARIABLES
    //----------------------------------------------------------------


    int maxMotifLength = 0;
    int numOfFactors = -1;
    ifstream transfac1(argv[1]);
    int m = 0;
    string row;
    if(!transfac1){
        cerr << "Matrix file not opened\n";
        exit(1);
    }

    while(!transfac1.eof()){
        getline(transfac1,row);
        if((row.substr(0,1) == "#")||(row.substr(0,1) == "")){
            continue;
        }
        if(row.substr(0,1) == ">"){
            numOfFactors++;
            if(m > maxMotifLength){
                maxMotifLength = m;
            }
            m = 0;
            continue;
        }
        m++;
    }

    if(m > maxMotifLength){
        maxMotifLength = m;
    }
    transfac1.close();

    maxMotifLength++;
    numOfFactors++;

    double lnR0[numOfFactors];
    string factornames[numOfFactors];
    int motiflength[numOfFactors];

    double *** complement;
    double *** matrix;
    matrix = new double ** [numOfFactors];
    complement = new double ** [numOfFactors];
    for(int i = 0; i < numOfFactors; i++){
        matrix[i] = new double * [maxMotifLength];
        complement[i] = new double * [maxMotifLength];
        for(int j = 0; j < maxMotifLength; j++){
            matrix[i][j] = new double[4];
            complement[i][j] = new double[4];
        }
    }


    //----------------------------------------------------------------
    //READ POSITION WEIGHT MATRICIES
    //----------------------------------------------------------------

    //reading file variables
    int factors = -1;
    string word[100]; //elements in row
    double max = 0; //consensus base count
    string delimiters = " \t"; //word seperators in each line
    int reading;
    int start, end;

    ifstream transfac(argv[1]);

    if(!transfac){
        cerr << "Matrix file not opened\n";
        exit(1);
    }

    while(!transfac.eof()){
        getline(transfac,row);
        start = row.find_first_not_of(delimiters);
        if(row.substr(0,1) == "#"){ // commentary line
            continue;
        }
        if(row.substr(0,1) == ""){ // empty line
            continue;
        }
        int i = 0;
        while(start != string::npos){ //split row into tokens - word[]
            end = row.find_first_of(delimiters,start + 1);
            if(end == string::npos){
                end = row.length();
            }
            word[i] = row.substr(start,end - start);
            i++;
            start = row.find_first_not_of(delimiters,end + 1);
        }

        if(i == 0){ // line without content
            continue;
        }

        if(word[0].substr(0,1) == ">"){ //new matrix reached
            factors++;
            motiflength[factors] = 0;
            factornames[factors] = word[1];//word[0].substr(2);//1
            for(int w = 0; w < i; w++){
                if(word[w] == "lnR0:"){
                    lnR0[factors] = strtod(word[w + 1].c_str(),NULL);
                }
            }
            continue;
        }
        matrix[factors][motiflength[factors]][A] = strtod(word[0].c_str(),NULL);
        matrix[factors][motiflength[factors]][C] = strtod(word[1].c_str(),NULL);
        matrix[factors][motiflength[factors]][G] = strtod(word[2].c_str(),NULL);
        matrix[factors][motiflength[factors]][T] = strtod(word[3].c_str(),NULL);
        motiflength[factors]++;
    }
    transfac.close();

    for(int f = 0; f <= factors; f++){
        for(int p = 0; p < motiflength[f]; p++){
            complement[f][motiflength[f]-p-1][A] = matrix[f][p][T];
            complement[f][motiflength[f]-p-1][T] = matrix[f][p][A];
            complement[f][motiflength[f]-p-1][G] = matrix[f][p][C];
            complement[f][motiflength[f]-p-1][C] = matrix[f][p][G];
        }
    }

    //----------------------------------------------------------------
    //READ FASTA FILE
    //----------------------------------------------------------------
    for(int i = 0; i <= factors; i++){
        cout << "\t" << factornames[i];
    }
    cout << "\n";
    ifstream fasta(argv[2]);

    if(!fasta){
        cerr << "FASTA file not opened\n";
        exit(1);
    }

    int last = 1; //EXIT CONDITION FOR LAST SEQUENCE IN FASTA FILE
    string newID, seqID;
    int seqlength = 0;
    int sequenceheader = 0;
    int firstrun = 1; //indicates first header
    string bases; //complete sequence
    string newbases; //sequence in each row of fasta file

    //OPEN OUTPUT SEQUENCES FILES
    ofstream file_output_sequences;
    file_output_sequences.open(argv[4]);
    file_output_sequences << "#TF transcription factor\n";
    file_output_sequences << "#AFFINITY_VALUE best affinity value predicted in SEQUENCE_NUMB (see filtered_sequences.fa for that)\n";
    file_output_sequences << "#LINE - 1-based index in filtered_sequences.fa\n";
    file_output_sequences << "#START_POS - 0-based start position in sequence\n";
    file_output_sequences << "#LENGTH - length of sequence (END_POS = START_POS + LENGTH)\n";
    file_output_sequences << "#CURRENT_AREA - current area on the actual chromosome\n";
    file_output_sequences << "TF\tAFFINITY_VALUE\tLINE\tSTART_POS\tLENGTH\tSEQUENCE\tCURRENT_AREA\n";

    //TODO: uncomment!!
    omp_set_num_threads(num_workers);
    int line_fasta = -1;
    std::string current_area ="";
    while(!fasta.eof())
    { //GO THROUGH FASTA FILE
        line_fasta++;
        getline(fasta,newbases);
        if(newbases.substr(0,1) == "#"){ //
            continue;
        }
        if(newbases.substr(0,1) == ">"){ //NEW SEQUENCE HEADER
            current_area=newbases;
            string delimiters = " \t"; //word seperators in each line
            int i = 0;
            start = newbases.find_first_not_of(delimiters);
            while(start != string::npos){ //split row into tokens - word[]
                end = newbases.find_first_of(delimiters,start + 1);
                if(end == string::npos){
                    end = newbases.length();
                }
                word[i] = newbases.substr(start,end - start);
                i++;
                start = newbases.find_first_not_of(delimiters,end + 1);
            }
            sequenceheader = 1;
            newID = word[0].substr(1);
        }

        if(sequenceheader == 0){ //SEQUENCE LINE
            bases += newbases;
            seqlength = bases.length();
            continue;
        }

        if(firstrun == 1){ //SKIP ANNOTATION FOR FIRST HEADER
            firstrun = 0;
            sequenceheader = 0;
            seqID = newID;
            continue;
        }

        label: int jumpLastSeq;
        cout << seqID ;
        //New: initialize array for all affinity values of one TF
        double affinityValuesTF[factors+1];
        int bestStartTF[factors+1];
        int bestLengthTF[factors+1];
        std::vector<std::string> bestSequenceTF(factors+1);
        std::vector<std::string> bestSequenceTF_name(factors+1);
        double bestStartTF_affinityValueTF[factors+1];
        //LOOP OVER FACTORS

        //TODO: UNCOMMENT!!
#pragma omp parallel for
        for(int f = 0; f <= factors; f++){
            affinityValuesTF[f] = 0;
            bestStartTF[f] = 0;
            bestLengthTF[f] = 0;
            bestStartTF_affinityValueTF[f] = DBL_MIN;
            bestSequenceTF_name[f] = factornames[f];

            int illegalBase, BASE;
            double P_combined = 0 ; //only palindrome correction, for entire seq
            double product, P_bound_F, P_bound_C,P_toBeAdded;
            double dE_forward, dE_compl;
            for(int n = 0; n < seqlength - motiflength[f] + 1; n++){ //LOOP OVER SEQUENCE
                dE_forward = 0;
                dE_compl = 0;

                std::string current_sequence = "";

                for(int m = 0; m < motiflength[f]; m++){ //LOOP OVER MOTIF
                    illegalBase = 0;
                    switch(bases[n + m])
                    {
                        case 'A':
                            BASE = 0;
                            current_sequence += "A";
                            break;
                        case 'C':
                            BASE = 1;
                            current_sequence += "C";
                            break;
                        case 'G':
                            BASE = 2;
                            current_sequence += "G";
                            break;
                        case 'T':
                            BASE = 3;
                            current_sequence += "T";
                            break;
                        case 'a':
                            BASE = 0;
                            current_sequence += "A";
                            break;
                        case 'c':
                            BASE = 1;
                            current_sequence += "C";
                            break;
                        case 'g':
                            BASE = 2;
                            current_sequence += "G";
                            break;
                        case 't':
                            BASE = 3;
                            current_sequence += "T";
                            break;
                        default:
                            illegalBase = 1;
                    }
                    if(illegalBase == 1){break;}
                    dE_forward += matrix[f][m][BASE];
                    dE_compl += complement[f][m][BASE];
                }//loop over motif

                //CALCULATE P(BOUND) FOR CURRENT SITE
                if(illegalBase == 0){
                    product = exp(lnR0[f] - dE_forward);
                    P_bound_F = product/(1 + product);
                    product = exp(lnR0[f] - dE_compl);
                    P_bound_C = product/(1 + product);
                    P_toBeAdded = P_bound_F + (1 - P_bound_F) * P_bound_C;

                    //check if sequence has higher affinity than other sequences in this area
                    if (P_toBeAdded > bestStartTF_affinityValueTF[f])
                    {
                        bestStartTF[f] = n;
                        bestSequenceTF.at(f) = current_sequence;
                        bestLengthTF[f] = motiflength[f];
                        bestStartTF_affinityValueTF[f] = P_toBeAdded;
                    }

                    P_combined += P_bound_F + (1 - P_bound_F) * P_bound_C;
                }
            }//loop over sequence
            affinityValuesTF[f]= P_combined;
        }//loop over factors

        //Print the affinity values for the current sequence of all TFs
        for(int f = 0; f <= factors; f++){
            cout << "\t" << setprecision(10) << affinityValuesTF[f];

            //file_output_sequences << "TF\tAFFINITY_VALUE\tSEQUENCE_NUMB\tSTART_POS\tLENGTH\tSEQUENCE\tCURRENT_AREA\n";
            file_output_sequences << bestSequenceTF_name[f] << "\t" << affinityValuesTF[f] << "\t" << line_fasta << "\t" << bestStartTF[f] << "\t" << bestLengthTF[f] << "\t" << bestSequenceTF.at(f) << "\t" << current_area <<"\n";

        }
        cout << "\n";

        //RESET SEQUENCE VARIABLES
        sequenceheader = 0;
        seqID = newID;
        bases = "";
    }//loop over fasta file

    fasta.close();


    //LAST SEQUENCE
    if(last == 1){
        last = 0;
        goto label;
    }

    //CLOSE OUTPUT SEQUENCES FILES
    file_output_sequences.close();

    for(int i = 0; i < numOfFactors; i++){
        for(int j = 0; j < maxMotifLength; j++){
            delete [] matrix[i][j];
            delete [] complement[i][j];
        }
        delete [] matrix[i];
        delete [] complement[i];
    }
    delete [] matrix;
    delete [] complement;
    return 0;
}
