#include <unordered_map>
#include <iostream>
#include <vector>
#include <queue>
#include <bitset>

using namespace std;

/*
**************************
1. Build a Frequency Table
**************************
- i.e. Frequency of each pixel value in the image. 
    (for grayscale images, we count of how often each pixel intensity (0-255) appears) 

Steps for Frequency Table:
    - Loop over all the pixels in the image.
    - Count how often each value (or group of values) appears.

Details:
    - assumes image = 1D vector (buffer), where each pixel = 8-bit value
    - The unordered_map freqTable maps pixel value to its frequency in the image
*/


void buildFrequencyTable(const vector<uint8_t>& buffer, unordered_map<uint8_t, int>& freqTable) {
    // Traverse the image buffer and count the frequency of each pixel value
    for (size_t i = 0; i < buffer.size(); i++) {
        freqTable[buffer[i]]++;  // Increment the frequency of the pixel value
    }
}


/*
*************************
2. Generate Huffman Codes
*************************
    - frequency table is used to build a Huffman tree. 
    - Tree built by repeatedly combining the least frequent pixel values into larger nodes until we have one tree

Steps:
    - Create a Priority Queue (i.e. Min-Heap): efficiently selects least frequent elements
    - Build the Huffman Tree: Combine nodes with lowest frequencies until only one tree remains
    - Assign Codes: Traverse Huffman tree and assign binary codes based on tree structure

Details:
    - Build Huffman tree using a priority queue (std::priority_queue), starting with the least frequent nodes
    - The generateHuffmanCodes() function traverses Huffman tree and generates binary codes for each pixel, and stores codes in huffCodes map.
*/

// Huffman tree node
struct HuffmanNode {
    uint8_t value;  // Pixel value (for grayscale)
    int freq;       // Frequency of this pixel value
    HuffmanNode* left;
    HuffmanNode* right;

    HuffmanNode(uint8_t v, int f) : value(v), freq(f), left(nullptr), right(nullptr) {}
};

// Comparison function to build a priority queue (min-heap)
struct Compare {
    bool operator()(HuffmanNode* left, HuffmanNode* right) {
        return left->freq > right->freq;
    }
};

// Build the Huffman tree
HuffmanNode* buildHuffmanTree(const unordered_map<uint8_t, int>& freqTable) {
    priority_queue<HuffmanNode*, vector<HuffmanNode*>, Compare> pq;

    // Create a leaf node for each pixel value and add to the priority queue
    for (const auto& pair : freqTable) {
        pq.push(new HuffmanNode(pair.first, pair.second));
    }

    // Build the tree by combining the two nodes with the lowest frequency
    while (pq.size() > 1) {
        HuffmanNode* left = pq.top(); pq.pop();
        HuffmanNode* right = pq.top(); pq.pop();

        HuffmanNode* parent = new HuffmanNode(0, left->freq + right->freq);
        parent->left = left;
        parent->right = right;

        pq.push(parent);
    }

    // The root of the tree
    return pq.top();
}

// Generate Huffman codes from the tree
void generateHuffmanCodes(HuffmanNode* root, string str, unordered_map<uint8_t, string>& huffCodes) {
    if (!root) return;

    // If we reached a leaf node, assign the Huffman code
    if (!root->left && !root->right) {
        huffCodes[root->value] = str;
    }

    // Traverse left and right subtrees
    generateHuffmanCodes(root->left, str + "0", huffCodes);
    generateHuffmanCodes(root->right, str + "1", huffCodes);
}


/*
*******************
3. Encode the Image
*******************
Replace pixel values with their corresponding Huffman codes to create bitstream

Details:
    - Converts the image’s pixel data (buffer) into single bitstream based on Huffman codes stored in "huffCodes"
*/

string encodeImage(const vector<uint8_t>& buffer, const unordered_map<uint8_t, string>& huffCodes) {
    string encodedStr;
    
    // For each pixel, append its Huffman code to the encoded string
    for (size_t i = 0; i < buffer.size(); i++) {
        encodedStr += huffCodes.at(buffer[i]);
    }

    return encodedStr;
}

/*
****************************
4. Write the Compressed Data
****************************

Writing compressed Huffman data to TIFF file involves modifying the TIFF structure to handle raw compressed bitstreams, which is not straightforward and sometimes difficult. 
Additionally, most standard TIFF libraries (e.g. libtiff) don't directly support raw Huffman-encoded data. 

Therefore, we need to:
    - Write custom TIFF header file
    - Store the compressed bitstream in some form of a TIFF data block
    - Ensure decompression logic (Huffman decoding) is available when reading the file back
    - possibly store the compressed data in an external file (like a .bin or .huff file), and later decode it when needed

A hypothetical approach to store compressed data in TIFF:

    - modify TIFF writer to store the Huffman-encoded data as raw bytes
    - custom compression type in TIFF header (i.e. TIFFTAG_COMPRESSION), and define it as "Huffman" or use "custom" compression scheme

Summary of Steps:
    - Build frequency table of pixel values
    - Generate Huffman tree using a frequency table
    - Encode image by replacing pixel values with corresponding Huffman codes
    - Write compressed data into a TIFF file, modifying TIFF format to support raw Huffman data, and/or using a simpler external file storage
Conclusion:
While expanding on our huffman encoding homework by implementing Huffman compression is a good learning exercise, the realities of integrating this into the TIFF format complex, and beyond the scope of this project timeline.  Furthermore, from my research, a TIFF is not designed for raw bitstream compression like Huffman. Additional expansion of this project scope could involve including libraries like zlib (utilizing the DEFLATE feature, which uses Huffman coding), or just using JPEG compression for practical use cases of analysis of my code's ability to work with several types of images.
*/
