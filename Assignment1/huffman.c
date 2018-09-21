/* ---------------------------------------
 * Description: COMP9319 Assignment 1
 * Author: Changfeng LI(Frank)(z5137858)
 * file: huffman.c
 * Recently Updated: 2018-08-24 19:10
 * Version: 2.0
 * ---------------------------------------*/

/* ==========
// current progress:
!!!perfectly done encoding
!!!perfectly done decoding
!!!done search
=== SUBMIT
=========== */

/*--------------------------------------------------------*/
/*                    Headfile including                  */
/*--------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
//#include <malloc.h>
//#include "windows.h"


/*--------------------------------------------------------*/
/*                    Macro Definition                    */
/*--------------------------------------------------------*/
#define MAX_HEADER_BYTE 1024
#define SIZE_ALPHABET 256 // size of all possible char
// worst case of len of huffman code when tree grow to one side
// the len = nb_feat = 256
#define MAX_TREE_HT SIZE_ALPHABET
#define EVENT_FAILURE 1
#define EVENT_SUCCESS 0
#define END_OF_FILE -1
#define NB_CHAR_BIT 8
#define SIZE_BIT_BUFFER NB_CHAR_BIT
#define BIT_VOID 2
#define BIT_TRUE 1
#define BIT_FLSE 0

#define GET_ARRAY_LEN(arr, len) {len = (sizeof(arr) / sizeof(arr[0]));}


/*--------------------------------------------------------*/
/*                    Type Definition                     */
/*--------------------------------------------------------*/
//B_Tree node
typedef struct nodet {
    short index;
    unsigned int weight;
    struct nodet *lchild, *rchild;
} NodeT, *BTree;

//Hash_Table node
typedef struct nodeh {
    unsigned char sKey;
    unsigned int* codeArr;
    unsigned int codeLen;
} NodeH;

//Header which saves information to exactly 1024 Bytes
typedef struct {
    //int len; // 4 Bytes
    unsigned int freqArr[SIZE_ALPHABET]; // 1024 Bytes
    // unsigned short bitArr[MAX_HEADER_BYTE]; // 1024 Bytes
    // char preOrder[MAX_HEADER_BYTE]; // 1024 Bytes
    // char midOrder[MAX_HEADER_BYTE]; // 1024 Bytes
    // char sufOrder[MAX_HEADER_BYTE]; // 1024 Bytes
    // char poem[MAX_HEADER_BYTE]; //1024 Bytes
} Header;


/*--------------------------------------------------------*/
/*                    Global Variable                     */
/*--------------------------------------------------------*/
unsigned int nb_feat_appear = 0;
// max size of test file is 100MB, if it contains single char,
// then freq will be 100_000_000 / 1 = 100_000_000, so use uint
unsigned int freq[SIZE_ALPHABET];
NodeT *nodes = NULL;
unsigned int nb_nodes = 0;
int cursor = 1;
char *s = "";
unsigned char *buffer;
unsigned int bit_table[9] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0xFF};
int bits_in_buffer = 0;
int current_bit = 0;
int flag_eof = 0;
unsigned int original_size = 0;
int vec_left[MAX_TREE_HT] = {0};
unsigned int nb_iter = 0;
unsigned int *g ;
char ch_global;
NodeH* hashTable[SIZE_ALPHABET]; //hash table data strcutrue
int hash_table_size;  //the number of key-value pairs in the hash table!


/*--------------------------------------------------------*/
/*                 Function Declaration                   */
/*--------------------------------------------------------*/
short fileVerify(const char* i_f, const char* o_f, unsigned char tag) ;
void getFreqTable(FILE *fp);
void generateTree(void);
void createNode(int index, int weight, int index_lchild, int index_rchild);
int readInputFile(FILE *f);
int writeOutputFile(FILE *fi, FILE *fo, char file_class);
int readBit(FILE *f, unsigned char finish_tag);
int writeBit(FILE *f, unsigned int bit, unsigned char finish_tag);
void printArrayNodeT(NodeT *array);
void printArrayInt(int *array, int len);
void displayTree(NodeT* btroot, int identity);
unsigned int hash_table_hash_str(const char skey);
void hashTableInsert(const unsigned char skey, unsigned int* codearr, unsigned int codelen);
NodeH* hashTableFindItem(const unsigned char skey);
void hashFree();
void printHash();
short huffman_encoder(const char* i_f, const char *o_f);
void huffmanEncoding(NodeT *btroot);
void getCodes(NodeT *btroot, unsigned int arr[],  char bit_arr[], int depth);
short huffman_decoder(const char *i_f, const char *o_f);
void huffmanDecoding(FILE *fi, FILE *fo, char file_class);
short huffman_searcher(char *query, const char *i_f);
unsigned int huffmanSearching(char *q);
char *append(const char *oldstring, const char c);
unsigned int matchString(char *substr, char *oristr);
void dropSys();
void initSys();

/*--------------------------------------------------------*/
/*                  Essential Function                    */
/*--------------------------------------------------------*/
// ==================
// Function: init
//
// ==================
void initSys() {
    // nodes array init -- only allocate for maximum possible child nodes
    nodes = (NodeT *) calloc(SIZE_ALPHABET, sizeof(NodeT));
    // hash table init
    hash_table_size = 0;
    memset(hashTable, 0, sizeof(NodeH *) * SIZE_ALPHABET);
    buffer = (unsigned char *) calloc (1, sizeof(unsigned char));

}

// ==================
// Function: drop
//
// ==================
void dropSys() {
    free(nodes);
    free(buffer);
    hashFree();
}

short fileVerify( const char* i_f, const char* o_f, unsigned char tag) {
    FILE *f_i, *f_o;
    if (tag == 'e' || tag == 'd') {
        if ((f_i = fopen(i_f, "rb")) == NULL) {
            perror("Failed to open input file");
            fclose(f_i);
            return EVENT_FAILURE;
        }
        if ((f_o = fopen(o_f, "wb")) == NULL) {
            perror("Failed to open output file");
            fclose(f_o);
            return EVENT_FAILURE;
        }
    }
    else if (tag == 's'){
        if ((f_i = fopen(i_f, "rb")) == NULL) {
            perror("Failed to open input file");
            fclose(f_i);
            return EVENT_FAILURE;
        }
    }
    else {
        return EVENT_FAILURE;
    }
    return EVENT_SUCCESS;
}
// ==================
// Function: getFreqTable
//
// ==================
void getFreqTable(FILE *fp) {
    int ch, tmp;

    while ((ch = fgetc(fp)) != EOF) {
        freq[ch]+=1;
    }
    for (tmp = SIZE_ALPHABET-1; tmp>=0; tmp--) {
        if (freq[tmp] > 0) {
            //printf("ch: %c ||| ascii: %d ||| freq: %d\n",ch,ch,freq[ch]);
            nb_feat_appear += 1;
        }
    }
}

// ==================
// Function: generateTree
//
// ==================
void generateTree() {
    int i, curr_freq, tag, il, ir;

    // Step1: realloc space
    // 2n - 1 many nodes overall + node at index 0 is always a empty node = 2n, this is useful for computing
    nodes = (NodeT *) calloc(2 * nb_feat_appear, sizeof(NodeT));
    // Step2: add all child nodes
    for (i = 0; i < SIZE_ALPHABET; ++i) {
        curr_freq = freq[i];
        if (curr_freq > 0) {
            ////printf("curr alpha = %c || curr index = %d || curr weight = %d\n", i, tag, curr_freq);
            tag = i;
            createNode(tag, curr_freq, -1, -1);
        }
    }

    // Step3: add parent nodes to build tree

    while (cursor < nb_nodes) {
        il = cursor;
	cursor++;
        ir = cursor;
	cursor++;
        curr_freq = nodes[il].weight + nodes[ir].weight;
        tag = -(ir/2);
        createNode(tag, curr_freq, il, ir);

    }
    //printArray(nodes);
    //displayTree(&nodes[num_nodes], 0);

}

void createNode(int index, int weight, int index_lchild, int index_rchild) {
    int j = nb_nodes++;
    while (weight < nodes[j].weight && j>0) {
	// ensure its ordering
        memcpy(&nodes[j + 1], &nodes[j], sizeof(NodeT));
        j-=1;
    }

    j+=1;
    nodes[j].weight = weight;
    nodes[j].index = index;
    // I assign !- as leaf nodes
    // I assign - as parent nodes
    if (nodes[j].index >= 0) {
        nodes[j].lchild = NULL;
        nodes[j].rchild = NULL;
    }
    else {
        nodes[j].lchild = &nodes[index_lchild];
        nodes[j].rchild = &nodes[index_rchild];
    }
}

// ================================
// ================================
// ================================
// ================================

// ******************* encoding
short huffman_encoder(const char* i_f, const char* o_f) {
    FILE *fi, *fo;

    if (fileVerify(i_f, o_f, 'e') == EVENT_FAILURE) {
        return EVENT_FAILURE;
    }
    else {
        fi = fopen(i_f, "rb");
        fo = fopen(o_f, "wb");
    }

    //printf("\n1.compute freq------\n\n");
    getFreqTable(fi);

    //printf("\n2.bulid_bitree------\n\n");
    generateTree();

    //printf("\n3.huff_encodes------\n\n");
    if (nb_nodes > 0)
    	huffmanEncoding(&nodes[nb_nodes]);

    //printf("\n4.write_file------\n\n");
    writeOutputFile(fi, fo, 'e');

    //unsigned long ori_size = ftell(fi);
    //unsigned long enc_size = ftell(fo);

    ////printf("==== before tar: size = %lu\n",ori_size);
    ////printf("==== after tar: size = %lu\n",enc_size);
    ////printf("==== percentage: size = %.8f\n",(float)(enc_size) / ori_size);

    fclose(fi);
    fclose(fo);
    free(g);

    return EVENT_SUCCESS;
}


void huffmanEncoding(NodeT *btroot)
{
    // get code and encoding huffman code to character by using alphabets
    unsigned int top = 0;
    unsigned int code_arr[MAX_TREE_HT];
    char bit_arr[MAX_HEADER_BYTE];

    getCodes(btroot, code_arr, bit_arr, top);

    bit_arr[nb_iter] = '\0';

    //printHash(); // we can print overall mid-order sequence if needed
    // for(int i=0;i<1024;i++) {
    //     //printf("%c,", bit_arr[i]);
    //     if (bit_arr[i] == '\0') {
    //         break;
    //     }

    // }
    // //printf("\n");
}


void getCodes(NodeT *btroot, unsigned int code_arr[], char bit_arr[], int depth) {
    // mid-order traversing the btree
    if (btroot->lchild)
    {
        code_arr[depth] = 1;//0
        bit_arr[nb_iter] = '1';
        nb_iter++;
        getCodes(btroot->lchild, code_arr, bit_arr, depth + 1);
    }

    if (btroot->rchild)
    {
        code_arr[depth] = 0;//1
        bit_arr[nb_iter] = '0';
        nb_iter++;
        getCodes(btroot->rchild, code_arr, bit_arr, depth + 1);
    }

    if (!btroot->lchild && !btroot->rchild)
    {

        g = (unsigned int *) calloc(depth, sizeof(unsigned int));

        //g[0] = depth;
	int t = 0;
        for(t=0; t<depth; t++) {
            g[t] = code_arr[t];

        }
        bit_arr[nb_iter] = (unsigned char)(btroot->index);
        nb_iter++;

        hashTableInsert((unsigned char)(btroot->index), g, depth);

    }
}


// ******************* decoding
short huffman_decoder(const char *i_f, const char *o_f) {

    FILE *fi, *fo;
    if (fileVerify(i_f, o_f, 'd') == EVENT_FAILURE) {
        return EVENT_FAILURE;
    }
    else {
        fi = fopen(i_f, "rb");
        fo = fopen(o_f, "wb");
    }
    //printf("\n1.read freq------\n\n");
    readInputFile(fi);

    //printf("\n2.bulid bitree------\n\n");
    generateTree();

    //printf("\n3.decode and gen original file------\n\n");
    huffmanDecoding(fi, fo, 'o'); //writeOutputFile(fi, fo, 'o');

    fclose(fi);
    fclose(fo);

    return EVENT_SUCCESS;
}

void huffmanDecoding(FILE *fi, FILE *fo, char file_class) {
    int i = 0;
    int bit;
    int node_index = nodes[nb_nodes].index;
    int root_index = node_index; // the last node which means it is a root node
    fseek(fi,0,SEEK_END);
    int file_size = ftell(fi);

    fseek(fi, MAX_HEADER_BYTE, SEEK_SET); // move to the head of input file.
    ////printf("===== original_size = %d\n",original_size);
    while (1) {
        bit = readBit(fi, (file_size==MAX_HEADER_BYTE)?'n':'r');

        if (bit == -1)//eof
            break;

        node_index = nodes[-(bit + node_index * 2)].index; // get the 2*index+(bit==1)?1:0 to get its l or r

        if (node_index >= 0) {
            char c = node_index;
	        if (file_class == 'o') {
		        ch_global = c;
		        writeOutputFile(fi, fo, 'o');

	        }
	        else {
		        s = append(s, c);
            }
            if (++i == original_size)
                break;
            node_index = root_index;
        }
    }
}

// ******************* searching

short huffman_searcher(char *query, const char *i_f) {
    FILE *fi;
    if (fileVerify(i_f, NULL, 's') == EVENT_FAILURE) {
        return EVENT_FAILURE;
    }
    else {
        fi = fopen(i_f, "rb");
    }

    //printf("\n1.read freq------\n\n");
    readInputFile(fi);

    //printf("\n2.bulid bitree------\n\n");
    generateTree();

    //printf("\n3.decode bit stream------\n\n");
    huffmanDecoding(fi, fi, 's');

    //printf("\n4.searching\n");
    printf("%u\n",huffmanSearching(query));

    ////printf("%s\n", s);

    fclose(fi);
    return EVENT_SUCCESS;
}

unsigned int huffmanSearching(char *q) {
    unsigned int pattern_match = matchString(q, s);
    return pattern_match;
}

// Stitching string with char
char *append(const char *oldstring, const char c) {
    int cal;
    char *newstring;
    cal = asprintf(&newstring, "%s%c", oldstring, c);
    if (cal == -1)
        newstring = NULL;
    return newstring;
}

// match string
unsigned int matchString(char *substr, char *oristr)
{
    int count = 0;
    char *s_tmp;
    char *o_tmp;
    if (*substr == '\0') {
	   return 0;
    }

    while(*oristr != '\0') {
        s_tmp = substr;
	    o_tmp = oristr;
        while((*o_tmp!='\0')&&(*s_tmp!='\0')&&(*o_tmp == *s_tmp)) {
            s_tmp++;
            o_tmp++;
        }
        if(*s_tmp=='\0') {
            ++count; // fully cmp all the bits, so this is a match
        }
        oristr++;
    }
    return count;

}


//============ bit operation ==========//

int readInputFile(FILE *f) {
    unsigned int ascii_index = 0;
    //size_t bytes_read;
    unsigned int *buff;

    buff = (unsigned int*) calloc(SIZE_ALPHABET, sizeof(unsigned int));
    // read 1024 Bytes to get freqArr
    while(fread(&buff[ascii_index],sizeof(unsigned int),sizeof(unsigned char),f)==1) {
        if (ascii_index > 255) {
            break;
        }

        freq[ascii_index] = buff[ascii_index];

        if (freq[ascii_index]>0) {
            nb_feat_appear++;
            original_size+=freq[ascii_index];

        }
        ascii_index++;

    }
    //nb_nodes = nb_feat_appear;
    free(buff);
    return EVENT_SUCCESS;
}

int readBit(FILE *f, unsigned char finish_tag) {
    if (finish_tag == 'r') {
        if (current_bit == bits_in_buffer) {
            if (flag_eof)
                return END_OF_FILE;
            else {
                size_t bytes_read = fread(buffer, 1, 1, f);
                
                if (bytes_read < 1 || feof(f)) {
                    flag_eof = 1;
                }
                bits_in_buffer = 8;
                current_bit = 0;
            }
        }
        int bit = (bool)(buffer[0] & bit_table[current_bit]);

        ++current_bit;
        return bit;
    }

    else if (finish_tag == 'n') {//null file
        return END_OF_FILE;
    }
    else {
        return END_OF_FILE;
    }

}

//=========
int writeBit(FILE *f, unsigned int bit, unsigned char finish_tag) {

    if (finish_tag == 'w') {
        if (bits_in_buffer == 8) {
            size_t bytes_written = fwrite(buffer, 1, 1, f);//min write char once a time
            //printf("===== %zu ==\n",bytes_written);
            if (ferror(f) || bytes_written < 1) return EVENT_FAILURE;
            memset(buffer, 0, 1);
            bits_in_buffer = 0;
        }

        (bit==BIT_TRUE)
        ?(buffer[0] |= bit_table[bits_in_buffer % 8])
        :(buffer[0] &= bit_table[8]);
        ++bits_in_buffer;
        return EVENT_SUCCESS;

    }
    else if (finish_tag == 'r') {
        //printf("... %d\n", bits_in_buffer);
        size_t bytes_written = fwrite(buffer, 1, 1, f);
        //printf("===== %zu ==\n",bytes_written);
        if (ferror(f) || bytes_written < 1) return EVENT_FAILURE;
        bits_in_buffer = 0;

        return EVENT_SUCCESS;
    }
    else {
        return EVENT_FAILURE;
    }

}

// use command to test
// $ xxd -b out.huffman
int writeOutputFile(FILE *fi, FILE *fo, char file_class) {
    // Step1: write header[1024]
    if (file_class == 'e') { // write encoded file
        Header header;
        memset(header.freqArr, 0, SIZE_ALPHABET);
	    unsigned short i=0;
        for (i=0;i<SIZE_ALPHABET;i++) {
            header.freqArr[i] = freq[i];
        }
        fwrite(&header.freqArr, SIZE_ALPHABET, sizeof(unsigned int), fo);

        // Step2: write encoding array
    	if (nb_nodes > 0) {
    		fseek(fi, 0, SEEK_SET); // move to the start of input file
    		int curr_char;
    		while ((curr_char = fgetc(fi)) != EOF) {
    		    int l=0;
    		    NodeH *code = NULL;

    		    code = hashTableFindItem(curr_char);
    		    
    		    while (l++<code->codeLen) {

    			    writeBit(fo, code->codeArr[l-1], 'w'); //write bit
    		    }
            }
            writeBit(fo, BIT_VOID, 'r');
    	}

    }
    else if (file_class == 'o') { // recover original file mode
        fwrite(&ch_global, 1, 1, fo);
    }
    else {
        return EVENT_FAILURE;
    }
    return EVENT_SUCCESS;

}


//============ hash fundamental
void hashTableInsert(const unsigned char skey, unsigned int* codearr, unsigned int codelen) {
    ////printf("--- start insert---\n");

    if(hash_table_size >= SIZE_ALPHABET) {
        printf("run out of hash table memory!\n");
        return;
    }

    unsigned short pos = (unsigned short)skey;
    NodeH* pNN = (NodeH*)calloc(1, sizeof(NodeH));
    memset(pNN, 0, sizeof(NodeH));

    pNN->sKey = skey;
    pNN->codeArr = (unsigned int*)calloc(codelen, sizeof(unsigned int));
    pNN->codeArr = codearr;
    pNN->codeLen = codelen;

    hashTable[pos] = pNN;

    hash_table_size++;

}

NodeH* hashTableFindItem(const unsigned char skey) {
    unsigned short pos = (unsigned short)skey;
    NodeH* pHead = hashTable[pos];

    if(hashTable[pos]) {
        if(pHead->sKey == skey) {
            return pHead;
        }
    }
    return NULL;
}

void hashFree() {
    int pos;
    for(pos = 0; pos < SIZE_ALPHABET; ++pos) {
        if(hashTable[pos]) {
            NodeH* pHead = hashTable[pos];
            NodeH* pTemp = pHead;
            if(pTemp) {
                free(pTemp);
            }
        }
    }
}


//===============================
// useful helper function
/***************
 * function: print node array
 ***************/
void printArrayNodeT(NodeT *array) {
    int length = nb_nodes;
    printf("length = %d\n", length);
    printf("[");
    int j = 0;
    for(j = 0; j <= nb_nodes; j++) {
        printf("(%d, %d), ", array[j].index, array[j].weight);
    }
    printf("]\n");
}

/***************
 * function: print int array
 ***************/
void printArrayInt(int *array, int len) {
    //int length = sizeof(array);
    //printf("length = %d\n", length);
    printf("[");
    int j = 0;
    for(j = 0; j < len; j++) {
        printf("\033[33m%d, ", array[j]);
    }
    printf("\033[0m]\n");
}

/***************
 * function: print hashtable
 ***************/
void printHash() {
    printf("===========content of hash table=================\n");
    int i;
    for(i = 0; i < SIZE_ALPHABET; ++i)

        if(hashTable[i]) {
            printf("%d=>", i);
            printf("%c:", hashTable[i]->sKey);
	    int j=0;
            for(j=0;j<hashTable[i]->codeLen;j++) {
                printf("%d, ", hashTable[i]->codeArr[j]);
            }
            printf("\n");
        }
}

/***************
 * function: print btree
 ***************/
void displayTree(NodeT * btroot, int identity) {
    if(identity > 0) {
        int i =0;
	for(i = 0; i < identity - 1; ++i) {
            printf(vec_left[i] ? "│   " : "    ");
        }
        printf(vec_left[identity-1] ? "├── " : "└── ");
    }

    if(! btroot) {
        printf("(null)\n");
        return;
    }

    printf("%d\n", btroot->weight);
    if(!btroot->lchild && !btroot->rchild) {
        return;
    }
    vec_left[identity] = 1;
    displayTree(btroot->lchild, identity + 1);
    vec_left[identity] = 0;
    displayTree(btroot->rchild, identity + 1);
}


// ============================== //
// ============================== //
// ======== Main Function ======= //
// ============================== //
// ============================== //
int main(int argc, char *argv[]) {
    if (argc != 4) { // judge if there are 4 args
        printf("Invalid number of command line input!\n");
        exit(1);
    }
    char *mode = argv[1]; //default
    char *src_input = argv[2];
    char *src_ouput = argv[3];

    initSys();
    // encode mode
    if ( strcmp(mode,"-e") == 0 ) {
        //printf(" +=+=+=+=+= current mode is encode! +=+=+=+=+=\n\n");
        huffman_encoder(src_input, src_ouput);
        //====================//
        /*
        //printf("insert testing.........\n");
        const char key1 = 'a';
        const char key2 = 'c';
        const char key3 = 'm';
        const unsigned int val1[] = {1,0,0};
        const unsigned int val2[] = {1,1,0};
        const unsigned int val3[] = {0,1};

        hashTableInsert(key1, val1, 3);
        hashTableInsert(key2, val2, 3);
        hashTableInsert(key3, val3, 2);
        */
        //====================//

    }
    // decode mode
    else if ( strcmp(mode,"-d") == 0 ) {
        //printf(" +=+=+=+=+= current mode is decode! +=+=+=+=+=\n\n");
        huffman_decoder(src_input, src_ouput);
    }
    // search mode
    else if ( strcmp(mode,"-s") == 0) {
        //printf(" +=+=+=+=+= current mode is search! +=+=+=+=+=\n\n");
        huffman_searcher(src_input,src_ouput);
    }
    else {
        //printf(" +=+=+=+=+= current mode is default +=+=+=+=+=\n\n");
    }
    dropSys();
    return EVENT_SUCCESS;
}


