#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cfloat>  // for FLT_MAX and -FLT_MAX

//#define std::uint32_t u32; 
//#define std::uint16_t u16;
//#define std::uint8_t u8;

// ******for brevity ********
//typedef std::uint32_t u32; 
//typedef std::uint16_t u16;
//typedef std::uint8_t u8;
// typdef instead of "#define", so its avoiding potential conflictions with other macros.  
// With typdef, this issue does not occur because it is a language feature rather than a preprocessor directive.
using namespace std;
// **************************

// Custom deleter for TIFF* to ensure proper cleanup
void TIFFCloseDeleter(TIFF* tiff) {
    if (tiff) {
        TIFFClose(tiff);
    }
}

// Function to read the TIFF image into a buffer
void readTiffImage(const char* filename, vector<uint8_t>& buffer, uint32_t& width, uint32_t& height, uint16_t& samples_per_pixel) {
    //TIFF* tiff = TIFFOpen(filename, "r");
    // Using std::unique_ptr with custom deleter to automatically close the TIFF file
    std::unique_ptr<TIFF, decltype(&TIFFCloseDeleter)> tiff(TIFFOpen(filename, "r"), TIFFCloseDeleter);
    
    if (!tiff) {
        cerr << "Could not open TIFF file!" << endl;
        return;
    }

    if (!TIFFGetField(tiff.get(), TIFFTAG_IMAGEWIDTH, &width) ||
        !TIFFGetField(tiff.get(), TIFFTAG_IMAGELENGTH, &height) ||
        !TIFFGetField(tiff.get(), TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel)) {
            cerr << "Error reading TIFF fields!" << endl;
            TIFFClose(tiff.get());
            return;
    }

    // Ensure width and height are even for sub-band splitting (wavelet transform)
    if (width % 2 != 0) {
        width++;
    }
    if (height % 2 != 0) {
        height++;
    }

    // Calculate the padded width and height to be a multiple of 4
    uint32_t paddedWidth = (width + 3) / 4 * 4;  // Round up to next multiple of 4
    uint32_t paddedHeight = (height + 3) / 4 * 4;  // Round up to next multiple of 4

    // Resize the buffer to the padded size
    buffer.resize(paddedWidth * paddedHeight * samples_per_pixel);

    // Read the image scanlines into the buffer
    for (uint32_t row = 0; row < height; row++) {
        TIFFReadScanline(tiff.get(), &buffer[row * paddedWidth * samples_per_pixel], row);
    }
    
}

// Daubechies-4 low-pass and high-pass filter coefficients
const float h0 = (1 + sqrt(3)) / (4 * sqrt(2));
const float h1 = (3 + sqrt(3)) / (4 * sqrt(2));
const float h2 = (3 - sqrt(3)) / (4 * sqrt(2));
const float h3 = (1 - sqrt(3)) / (4 * sqrt(2));
const float g0 = -h3;
const float g1 = h2;
const float g2 = -h1;
const float g3 = h0;

// Coefficients for the inverse Daubechies-4 transform:
const float ih0 = h0;
const float ih1 = -h1;
const float ih2 = h2;
const float ih3 = -h3;

const float ig0 = g0;
const float ig1 = -g1;
const float ig2 = g2;
const float ig3 = -g3;

// Function to perform 1D Daubechies-4 transform on an array
void daubechies1D(vector<float>& data) {
    int N = data.size();
    vector<float> temp(N);

    for (int i = 0; i < N / 2; i++) {
        temp[i] = h0 * data[2*i] + h1 * data[2*i+1] + h2 * data[2*i+2] + h3 * data[2*i+3];
        temp[i + N/2] = g0 * data[2*i] + g1 * data[2*i+1] + g2 * data[2*i+2] + g3 * data[2*i+3];
    }
    for (int i = 0; i < N; i++) {
        data[i] = temp[i];
    }
}

// 2D Daubechies-4 Transform: Apply 1D transform on rows and columns
void daubechies2D(vector<vector<float>>& image) {
    int rows = image.size();
    int cols = image[0].size();

    // Apply the 1D Daubechies transform on each row
    for (int i = 0; i < rows; i++) {
        daubechies1D(image[i]);
    }
    // Apply the 1D Daubechies transform on each column
    for (int j = 0; j < cols; j++) {
        vector<float> column(rows);
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        daubechies1D(column);
        for (int i = 0; i < rows; i++) {
            image[i][j] = column[i];
        }
    }
}

void applyBitShiftingQuantization(vector<float>& data, double sigma) {
    for (float& value : data) {
        // Soft thresholding: set small values to zero based on sigma
        if (abs(value) < sigma) {
            value = 0;
        } else {
            // Apply bit-shifting by first converting to integer, shifting, and then back to float
            int quantizedValue = static_cast<int>(value) >> 2;  // just decided to shift by 2 bits (keep it simple)
            value = static_cast<float>(quantizedValue);
        }
    }
}


// Function to apply wavelet compression with bit-shifting quantization
void waveletTransform(vector<uint8_t>& buffer, uint32_t& width, uint32_t& height, double sigma, int bitShiftAmount) {
    // Calculate the padded width and height to be a multiple of 4
    uint32_t paddedWidth = (width + 3) / 4 * 4;  // Round up to next multiple of 4
    uint32_t paddedHeight = (height + 3) / 4 * 4;  // Round up to next multiple of 4

    // Resize the buffer to the padded size
    buffer.resize(paddedWidth * paddedHeight);

    // Optionally, fill the padding with zeros
    fill(buffer.begin() + width * height, buffer.end(), 0);  // Padding with zero
    
    //Create temp buffer
    vector<double> tempBuffer(width * height);

    // Convert the uint8_t buffer to double for the wavelet transform
    for (uint32_t i = 0; i < width * height; i++) {
        tempBuffer[i] = static_cast<double>(buffer[i]);
    }

    // Create a 2D image matrix
    vector<vector<float>> image(height, vector<float>(width));
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            image[row][col] = static_cast<float>(tempBuffer[row * width + col]);
        }
    }

    // Apply the 2D Daubechies transform
    daubechies2D(image);

    // Apply bit-shifting quantization
    for (uint32_t row = 0; row < height; row++) {
        applyBitShiftingQuantization(image[row], bitShiftAmount);
    }

    // Convert back to 1D buffer
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            tempBuffer[row * width + col] = image[row][col];
        }
    }

    // Normalize to 0-255 range for output (8-bit)
    double minVal = *min_element(tempBuffer.begin(), tempBuffer.end());
    double maxVal = *max_element(tempBuffer.begin(), tempBuffer.end());

    if (maxVal == minVal) {
    // Handle the case where min and max are the same (b/c would this cause blank grey output?)
    fill(buffer.begin(), buffer.end(), 128);  // Arbitrary mid-level grayscale value
    } else {
        for (uint32_t i = 0; i < width * height; i++) {
            double normalized = (tempBuffer[i] - minVal) / (maxVal - minVal) * 255.0;
            buffer[i] = static_cast<uint8_t>(std::clamp(normalized, 0.0, 255.0));
        }
    }
}

// Function to write a TIFF image from a buffer
void writeTiffImage(const char* filename, const std::vector<uint8_t>& buffer, uint32_t width, uint32_t height, uint16_t samples_per_pixel) {
    std::unique_ptr<TIFF, decltype(&TIFFCloseDeleter)> tiff(TIFFOpen(filename, "w"), TIFFCloseDeleter);
    
    if (!tiff) {
        std::cerr << "Could not open TIFF file for writing!" << std::endl;
        return;
    }

    TIFFSetField(tiff.get(), TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff.get(), TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff.get(), TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    TIFFSetField(tiff.get(), TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tiff.get(), TIFFTAG_COMPRESSION, COMPRESSION_LZW);

    // Write image scanlines
    for (uint32_t row = 0; row < height; ++row) {
        TIFFWriteScanline(tiff.get(), const_cast<uint8_t*>(&buffer[row * width]), row);
    }

    // TIFF file auto-closed, when it goes out of scope
}

void inverseDaubechies1D(vector<float>& data) {
    int N = data.size();
    vector<float> temp(N);

    for (int i = 0; i < N / 2; i++) {
        temp[2 * i] = ih0 * data[i] + ih1 * data[i + N / 2];
        temp[2 * i + 1] = ih2 * data[i] + ih3 * data[i + N / 2];
    }

    for (int i = 0; i < N; i++) {
        data[i] = temp[i];
    }
}

void inverseDaubechies2D(vector<vector<float>>& image) {
    int rows = image.size();
    int cols = image[0].size();

    // Apply the inverse 1D Daubechies transform on each row
    for (int i = 0; i < rows; i++) {
        inverseDaubechies1D(image[i]);
    }

    // Apply the inverse 1D Daubechies transform on each column
    for (int j = 0; j < cols; j++) {
        vector<float> column(rows);
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        inverseDaubechies1D(column);
        for (int i = 0; i < rows; i++) {
            image[i][j] = column[i];
        }
    }
}

void inverseWaveletTransform(vector<uint8_t>& buffer, uint32_t width, uint32_t height) {
    // Convert the uint8_t buffer to double for inverse transform
    vector<double> tempBuffer(width * height);

    for (uint32_t i = 0; i < width * height; i++) {
        tempBuffer[i] = static_cast<double>(buffer[i]);
    }

    // Create a 2D image matrix for inverse transform
    vector<vector<float>> image(height, vector<float>(width));
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            image[row][col] = static_cast<float>(tempBuffer[row * width + col]);
        }
    }

    // Apply the inverse 2D Daubechies transform
    inverseDaubechies2D(image);

    // Convert back to 1D buffer after the inverse transform
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            tempBuffer[row * width + col] = image[row][col];
        }
    }

    // Normalize to 0-255 range for output (8-bit)
    double minVal = *min_element(tempBuffer.begin(), tempBuffer.end());
    double maxVal = *max_element(tempBuffer.begin(), tempBuffer.end());

    for (uint32_t i = 0; i < width * height; i++) {
        double normalized = (tempBuffer[i] - minVal) / (maxVal - minVal) * 255.0;
        buffer[i] = static_cast<uint8_t>(std::clamp(normalized, 0.0, 255.0));
    }
}

// Function to write a 2D matrix (sub-band) to a TIFF image
void writeSubbandToTiff(const char* filename, const vector<vector<float>>& subband, uint32_t width, uint32_t height) {
    vector<uint8_t> buffer(width * height);

    // Find the min and max values in the subband
    float minVal = FLT_MAX;
    float maxVal = -FLT_MAX;
    
    for (const auto& row : subband) {
        minVal = min(minVal, *min_element(row.begin(), row.end()));
        maxVal = max(maxVal, *max_element(row.begin(), row.end()));
    }

    // Avoid division by zero in case the min and max values are the same
    if (maxVal == minVal) {
        fill(buffer.begin(), buffer.end(), 128);  // Assign a mid-level grayscale value
    } else {
        // Normalize and convert the sub-band to 8-bit values (0-255)
        for (uint32_t i = 0; i < height; i++) {
            for (uint32_t j = 0; j < width; j++) {
                float normalized = ((subband[i][j] - minVal) / (maxVal - minVal)) * 255.0f;
                buffer[i * width + j] = static_cast<uint8_t>(std::clamp(normalized, 0.0f, 255.0f));
            }
        }
    }

    // Write the 8-bit grayscale image to TIFF
    writeTiffImage(filename, buffer, width, height, 1); // Assuming single channel (grayscale)
}

// Function to apply the Daubechies-4 wavelet transform on a 2D image and return sub-bands
void daubechies2DWithSubbands(vector<vector<float>>& image, 
                              vector<vector<float>>& LL, 
                              vector<vector<float>>& LH, 
                              vector<vector<float>>& HL, 
                              vector<vector<float>>& HH) {
    
    
    int rows = image.size();
    int cols = image[0].size();

    if (rows % 2 != 0 || cols % 2 != 0) {
        cerr << "Image dimensions must be even for sub-band splitting!" << endl;
        return;
    }

    // Temporary buffers to hold transformed rows and columns
    vector<vector<float>> tempRows(rows, vector<float>(cols));
    vector<vector<float>> tempCols(rows, vector<float>(cols));

    // Apply the 1D Daubechies transform on each row
    for (int i = 0; i < rows; i++) {
        daubechies1D(image[i]);  // Apply the 1D transform to each row
    }

    // Apply the 1D Daubechies transform on each column
    for (int j = 0; j < cols; j++) {
        vector<float> column(rows);
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        daubechies1D(column);  // Apply the 1D transform to each column
        for (int i = 0; i < rows; i++) {
            tempCols[i][j] = column[i];
        }
    }

    // Now, split the transformed image into four sub-bands: LL, LH, HL, HH
    int halfRows = rows / 2;
    int halfCols = cols / 2;

    // Approximation (LL)
    LL.resize(halfRows, vector<float>(halfCols));
    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            LL[i][j] = tempCols[i * 2][j * 2];
        }
    }

    // Horizontal Detail (LH)
    LH.resize(halfRows, vector<float>(halfCols));
    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            LH[i][j] = tempCols[i * 2][j * 2 + 1];
        }
    }

    // Vertical Detail (HL)
    HL.resize(halfRows, vector<float>(halfCols));
    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            HL[i][j] = tempCols[i * 2 + 1][j * 2];
        }
    }

    // Diagonal Detail (HH)
    HH.resize(halfRows, vector<float>(halfCols));
    for (int i = 0; i < halfRows; i++) {
        for (int j = 0; j < halfCols; j++) {
            HH[i][j] = tempCols[i * 2 + 1][j * 2 + 1];
        }
    }
}

int main() {
    // Input and output file paths
    const char* inputFile = "cameraman.tif";
    const char* outputFile = "wavlet_compressed_debug_sigma42_compressed.tif";
    const char* outputReconstructedFile = "wavlet_compressed_debug_sigma42_reconstructed.tif";
    
    uint32_t width, height = 0; // directly assign, to avoid potential memory errors
    uint16_t samples_per_pixel = 1;  // Assuming grayscale image (single sample per pixel)
    double sigma = 42;  // Thresholding value for quantization
    int bitShiftAmount = 2;  // Example: shift by 2 bits (divide by 8)

    // Step 1: Read the TIFF image into a buffer
    vector<uint8_t> buffer;
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Apply wavelet compression with bit-shifting quantization
    waveletTransform(buffer, width, height, sigma, bitShiftAmount);

    // Step 3: Write the compressed image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    // Optional Step: Decompose the image into sub-bands (LL, LH, HL, HH)
    // Convert to 2D matrix for easier manipulation
    vector<vector<float>> image(height, vector<float>(width));
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            image[row][col] = static_cast<float>(buffer[row * width + col]);
        }
    }

    vector<vector<float>> LL, LH, HL, HH;
    daubechies2DWithSubbands(image, LL, LH, HL, HH);

    // Optional: Save the sub-bands to TIFF files
    writeSubbandToTiff("LL_subband.tif", LL, width / 2, height / 2);  // Approximation sub-band
    writeSubbandToTiff("LH_subband.tif", LH, width / 2, height / 2);  // Horizontal detail sub-band
    writeSubbandToTiff("HL_subband.tif", HL, width / 2, height / 2);  // Vertical detail sub-band
    writeSubbandToTiff("HH_subband.tif", HH, width / 2, height / 2);  // Diagonal detail sub-band

    // Step 4: Apply the inverse wavelet transform to reconstruct the image
    inverseWaveletTransform(buffer, width, height);

    // Step 5: Write the reconstructed image to a new TIFF file
    writeTiffImage(outputReconstructedFile, buffer, width, height, samples_per_pixel);

    cout << "Wavelet compression with bit-shifting quantization and reconstruction completed!" << endl;
    return 0;
}